## function to sort the SKUs after the number of sales
function SORTSALES(trans::SparseMatrixCSC{Bool,Int64}, sku_weight::Vector{<:Real})
    sales = zeros(Int64, size(trans, 2), 4)
    sales[:, 3] = round.(Int64, sum(trans; dims = 1) ./ transpose(sku_weight))
    sales[:, 4] .= 0
    for i in axes(trans, 2)
        sales[i, 1] = i
    end
    sales = sortslices(sales; dims = 1, by = x->x[3], rev = true, alg = QuickSort)
    for i in axes(trans, 2)
        sales[i, 2] = i
    end
    # 1 column: index of SKU in allocation matrix
    # 2 column: index of SKU in sorted sales
    # 3 column: number of total sold SKU
    # 4 column: indicator whether SKU is allocated
    return sales::Matrix{Int64}
end

## function to sort the SKU pairs after the highest number of coappearances,
## weighted by the average storage requirement: p_ij = q_ij / ((s_i + s_j) / 2)
function SORTPAIRS(Q::AbstractMatrix{<:Real}, sku_weight::Vector{<:Real})
    pair_i = Int64[]
    pair_j = Int64[]
    pair_v = Float64[]
    if Q isa SparseMatrixCSC
        rows_q = rowvals(Q)
        vals_q = nonzeros(Q)
        for col in 1:size(Q, 2)
            for idx in nzrange(Q, col)
                row = rows_q[idx]
                if row > col
                    push!(pair_i, row)
                    push!(pair_j, col)
                    push!(pair_v, vals_q[idx] / ((sku_weight[row] + sku_weight[col]) / 2))
                end
            end
        end
    else
        @inbounds for i in 2:size(Q, 1)
            weight_i = sku_weight[i]
            for j in 1:(i - 1)
                push!(pair_i, i)
                push!(pair_j, j)
                push!(pair_v, Q[i, j] / ((weight_i + sku_weight[j]) / 2))
            end
        end
    end
    order = sortperm(pair_v; rev = true)
    pairs = Matrix{Int64}(undef, length(order), 2)
    pairs[:, 1] = pair_i[order]
    pairs[:, 2] = pair_j[order]
    return pairs::Matrix{Int64}
end

## function to check the number of coappearances in the selected warehouse
## for the selected product
function CURRENTCOAPP(X::Matrix{Bool}, warehouse::Int64, sku::Int64, Q::Matrix{<:Real})
    sum = 0.0
    @inbounds @simd for i in axes(X, 1)
        sum += X[i, warehouse] * Q[sku, i]
    end
    return sum
end

function CURRENTCOAPP(
    X::Matrix{Bool}, warehouse::Int64, sku::Int64, Q::SparseMatrixCSC{<:Real}
)
    out = 0.0
    rows = rowvals(Q)
    vals = nonzeros(Q)
    @inbounds for idx in nzrange(Q, sku)
        out += X[rows[idx], warehouse] * vals[idx]
    end
    return out
end

## function to start place the seeds in the greedy seeds heuristic
function GREEDYSEEDSTART!(
    sales::Matrix{Int64},
    X::Matrix{Bool},
    Q::AbstractMatrix{<:Real},
    capacity_left::Vector{<:Real},
    sku_weight::Vector{<:Real},
)
    for k in 2:size(capacity_left, 1)
        if capacity_left[k] > 0
            # list all SKUs that have not been allocated yet and
            top_sales = sales[sales[:, 4] .== 0, :]
            # discard from the list those SKUs that are not in the top decile
            top_sales = top_sales[1:max(1, round(Int64, size(top_sales, 1) / 10)), :]
            # find the SKU from the list with the least coappearances to the
            # seeds already placed in the other warehouses
            coapp = zeros(Float64, size(top_sales, 1))
            for i in 1:size(top_sales, 1)
                sku = top_sales[i, 1]
                for j in 1:size(capacity_left, 1)
                    coapp[i] += CURRENTCOAPP(X, j, sku, Q)
                end
                coapp[i] /= sku_weight[sku]
            end
            # allocate the SKU to the DC
            seed_sku = top_sales[findmin(coapp)[2], 1]
            if capacity_left[k] < sku_weight[seed_sku]
                SWAPREPAIR!(X, capacity_left, sku_weight, seed_sku)
            else
                X[seed_sku, k] = 1
                capacity_left[k] -= sku_weight[seed_sku]
            end
            sales[top_sales[findmin(coapp)[2], 2], 4] = 1
        end
    end
end

## function to allocate the SKU with the highest potential
function GREEDYSEEDMAIN!(
    sales::Matrix{Int64},
    X::Matrix{Bool},
    Q::AbstractMatrix{<:Real},
    capacity_left::Vector{<:Real},
    sku_weight::Vector{<:Real},
)
    for i in axes(Q, 1)
        if sales[i, 4] == 0
            best_allocation = zeros(Float64, size(capacity_left, 1))
            sku = sales[i, 1]
            for d in axes(capacity_left, 1)
                if capacity_left[d] >= sku_weight[sku]
                    # retrieve the list of SKUs allocated to d so far
                    # calculate the average of coapperances with the SKUs in d
                    best_allocation[d] = CURRENTCOAPP(X, d, sku, Q)/sku_weight[sku]
                end
            end
            if findmax(best_allocation)[1] > 0
                X[sales[i, 1], findmax(best_allocation)[2]] = 1
                capacity_left[findmax(best_allocation)[2]] -= sku_weight[sales[i, 1]]
            else
                d_best = findmax(capacity_left)[2]
                if capacity_left[d_best] < sku_weight[sales[i, 1]]
                    SWAPREPAIR!(X, capacity_left, sku_weight, sales[i, 1])
                else
                    X[sales[i, 1], d_best] = 1
                    capacity_left[d_best] -= sku_weight[sales[i, 1]]
                end
            end
            sales[i, 4] = 1
        end
    end
end

function GREEDYPAIRSMAIN!(
    pairs::Matrix{Int64},
    X::Matrix{Bool},
    capacity_left::Vector{<:Real},
    sku_weight::Vector{<:Real},
)
    n_allocated = 0
    for i in axes(pairs, 1)
        ## For each pair of SKUs in sorted list of pairs do
        for j in 1:2
            for d in 1:size(capacity_left, 1)
                if capacity_left[d] >= sku_weight[pairs[i, j]]
                    if sum(X[pairs[i, j], :]) == 0
                        X[pairs[i, j], d] = 1
                        capacity_left[d] -= sku_weight[pairs[i, j]]
                        n_allocated += 1
                        break
                    end
                end
            end
        end
        if n_allocated == size(X, 1)
            break
        end
    end
end

# function to allocate the SKUs in the main part of the Greedy Orders Heuristic
function GREEDYORDERSMAIN!(
    order_sort::Vector{Int64},
    trans::SparseMatrixCSC{Bool,Int64},
    capacity_left::Vector{<:Real},
    X::Matrix{Bool},
    sku_weight::Vector{<:Real},
)
    for order in order_sort
        for sku in axes(trans, 2)
            # Determine the warehouse d with the lowest index that has one unallocated space
            if trans[order, sku] == 1
                for d in 1:length(capacity_left)
                    if capacity_left[d] >= sku_weight[sku]
                        if sum(X[sku, :]) == 0
                            X[sku, d] = 1
                            capacity_left[d] -= sku_weight[sku]
                        end
                    end
                end
            end
        end
    end
    for sku in axes(trans, 2)
        if sum(X[sku, :]) == 0
            d_best = argmax(capacity_left)
            if capacity_left[d_best] < sku_weight[sku]
                SWAPREPAIR!(X, capacity_left, sku_weight, sku)
            else
                X[sku, d_best] = 1
                capacity_left[d_best] -= sku_weight[sku]
            end
        end
    end
end

## function to calculate the number of SKUs that can be allocated
## multiple times
function BESTSELLING_B(capacity_left::Vector{<:Real}, sku_weights::Vector{<:Real})
    floor(Int64, (sum(capacity_left)-sum(sku_weights))/(size(capacity_left, 1)-1))
end

## function to initialise the bestselling heuristic
function BESTSELLINGSTART!(
    sales::Matrix{Int64},
    capacity_left::Vector{<:Real},
    B::Int64,
    X::Matrix{Bool},
    sku_weight::Vector{<:Real},
)
    nonuniform = !all(w -> w == sku_weight[1], sku_weight)
    if nonuniform
        remaining_w = sort(
            [sku_weight[sales[j, 1]] for j in axes(sales, 1) if sales[j, 4] == 0];
            rev = true,
        )
    end
    for d in axes(capacity_left, 1)
        if capacity_left[d] < B
            i = 1
            while i <= size(sales, 1) && capacity_left[d] >= sku_weight[sales[i, 1]]
                # the top sellers are replicated into every small warehouse;
                # only skip SKUs already placed in this particular warehouse
                if X[sales[i, 1], d] == 1
                    i += 1
                    continue
                end
                w_i = sku_weight[sales[i, 1]]
                # Check feasibility before committing (non-uniform weights only);
                # remaining_w tracks unallocated SKUs, so only remove a SKU on
                # its first placement, not when it is replicated
                if nonuniform
                    first_placement = sales[sales[i, 2], 4] == 0
                    capacity_left[d] -= w_i
                    idx_rm = 0
                    if first_placement
                        idx_rm = searchsortedfirst(remaining_w, w_i; rev = true)
                        deleteat!(remaining_w, idx_rm)
                    end
                    if !FEASIBLE_REMAINING(capacity_left, remaining_w)
                        capacity_left[d] += w_i
                        if first_placement
                            insert!(remaining_w, idx_rm, w_i)
                        end
                        break
                    end
                else
                    capacity_left[d] -= w_i
                end
                X[sales[i, 1], d] = 1
                sales[sales[i, 2], 4] = 1
                i += 1
            end
        end
        B = BESTSELLING_B(capacity_left, sku_weight)
    end
end

## function to allocate the SKUs that have been sold most in each
## warehouse until the remaining SKUs can only be allocated once
function BESTSELLINGTOP!(
    sales::Matrix{Int64},
    capacity_left::Vector{<:Real},
    capacity::Vector{<:Real},
    B::Int64,
    X::Matrix{Bool},
    sku_weight::Vector{<:Real},
)
    nonuniform = !all(w -> w == sku_weight[1], sku_weight)
    if nonuniform
        remaining_w = sort(
            [sku_weight[sales[j, 1]] for j in axes(sales, 1) if sales[j, 4] == 0];
            rev = true,
        )
    end
    for d in axes(capacity_left, 1)
        i = 1
        while capacity_left[d] > capacity[d] - B &&
                  i <= size(sales, 1) &&
                  capacity_left[d] >= sku_weight[sales[i, 1]]
            # the top B sellers are replicated into every warehouse;
            # only skip SKUs already placed in this particular warehouse
            if X[sales[i, 1], d] == 1
                i += 1
                continue
            end
            w_i = sku_weight[sales[i, 1]]
            # Check feasibility before committing (non-uniform weights only);
            # remaining_w tracks unallocated SKUs, so only remove a SKU on
            # its first placement, not when it is replicated
            if nonuniform
                first_placement = sales[sales[i, 2], 4] == 0
                capacity_left[d] -= w_i
                idx_rm = 0
                if first_placement
                    idx_rm = searchsortedfirst(remaining_w, w_i; rev = true)
                    deleteat!(remaining_w, idx_rm)
                end
                if !FEASIBLE_REMAINING(capacity_left, remaining_w)
                    capacity_left[d] += w_i
                    if first_placement
                        insert!(remaining_w, idx_rm, w_i)
                    end
                    break
                end
            else
                capacity_left[d] -= w_i
            end
            X[sales[i, 1], d] = 1
            sales[sales[i, 2], 4] = 1
            i += 1
        end
    end
end
