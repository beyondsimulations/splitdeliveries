# Check if remaining unallocated SKUs can be feasibly assigned to warehouses.
# Simulates greedy assignment (heaviest first, argmax warehouse).
# Used during replication phases to prevent over-replication.
function FEASIBLE_REMAINING(cap_left::Vector{<:Real}, remaining_weights::Vector{<:Real})
    sim_cap = copy(cap_left)
    for w in sort(remaining_weights; rev = true)
        d = argmax(sim_cap)
        sim_cap[d] < w && return false
        sim_cap[d] -= w
    end
    return true
end

# Local swap repair: when a heavy SKU can't fit anywhere, evict lightest
# single-assigned SKU(s) from a warehouse, place the heavy one, and
# re-allocate the evicted SKUs elsewhere. Tries each warehouse in order
# of decreasing capacity and only commits if evicted SKUs can be re-placed.
# Returns (d_placed, evicted_skus, d_old_list, d_new_list) for state repair.
function SWAPREPAIR!(
    X::Matrix{Bool}, cap_left::Vector{<:Real}, sku_weight::Vector{<:Real}, heavy_sku::Int64
)
    w_heavy = sku_weight[heavy_sku]

    for d in sortperm(cap_left; rev = true)
        # Find single-assigned SKUs in d (safe to move), sorted lightest first
        candidates = [
            (j, sku_weight[j]) for
            j in 1:size(X, 1) if X[j, d] && sum(@view(X[j, :])) == 1 && j != heavy_sku
        ]
        sort!(candidates; by = x->x[2])

        # Determine which SKUs to evict
        freed = 0.0
        to_evict = Int[]
        for (j, w) in candidates
            cap_left[d] + freed >= w_heavy && break
            freed += w
            push!(to_evict, j)
        end
        cap_left[d] + freed < w_heavy && continue

        # Simulate: can evicted SKUs be re-placed after heavy placement?
        sim_cap = copy(cap_left)
        sim_cap[d] += freed - w_heavy
        evicted_w = [sku_weight[j] for j in to_evict]
        FEASIBLE_REMAINING(sim_cap, evicted_w) || continue

        # Commit: evict, place heavy, re-allocate evicted
        for j in to_evict
            X[j, d] = false
            cap_left[d] += sku_weight[j]
        end
        X[heavy_sku, d] = true
        cap_left[d] -= w_heavy

        d_old = fill(d, length(to_evict))
        d_new = zeros(Int, length(to_evict))
        for (idx, j) in enumerate(to_evict)
            d2 = argmax(cap_left)
            X[j, d2] = true
            cap_left[d2] -= sku_weight[j]
            d_new[idx] = d2
        end
        return d, to_evict, d_old, d_new
    end

    error(
        "SWAPREPAIR!: Cannot fit SKU $heavy_sku (weight $w_heavy) in any warehouse via swap.",
    )
end

# Function to calculate the Coappearance Matrix Q
function COAPPEARENCE(trans::SparseMatrixCSC{Bool,Int64}, sku_weight::Vector{<:Real})
    Q64 = trans' * trans
    Q = SparseMatrixCSC{Int32,Int32}(
        size(Q64, 1),
        size(Q64, 2),
        Int32.(Q64.colptr),
        Int32.(rowvals(Q64)),
        Int32.(nonzeros(Q64)),
    )
    return Q
end

# function to set all entries on the principle diagonal to zero
function CLEANPRINCIPLE!(Q::AbstractMatrix{<:Real})
    for i in axes(Q, 1)
        Q[i, i] = 0
    end
    if Q isa SparseMatrixCSC
        dropzeros!(Q)
    end
end

# Functions to evaluate the number of parcels necessary to fullfil all orders
## COMBINEWAREHOUSES: create all possible warehouse combinations 
function COMBINEWAREHOUSES(capacity::Array{Int64,1})
    combination = [[k] for k in axes(capacity, 1)]
    combinations_X = combinations(combination)
    combination = [c for c in combinations_X]
    return combination::Array{Array{Array{Int64,1},1},1}
end

## PARCELSSEND: number of parcels due to ordersplitting necessary to fulfill all orders
function PARCELSSEND(
    trans::SparseMatrixCSC{Bool,Int64},
    X::Array{Bool,2},
    capacity::Array{Int64,1},
    combination::Array{Array{Array{Int64,1},1},1},
)
    if !any(iszero, sum(X; dims = 2))
        if sum(capacity) == size(trans, 2)
            parcel = trans * X
            parcel = sum(parcel .> 0)
        else
            warehouse_combination = zeros(Int64, size(X, 1), size(combination, 1))
            parcels_out = zeros(size(combination, 1))
            for i in axes(combination, 1)
                for j in axes(combination[i], 1)
                    for k in axes(X, 1)
                        if X[k, combination[i][j]] == [1]
                            warehouse_combination[k, i] = 1
                        end
                    end
                    parcels_out[i] += 1
                end
            end
            warehouse_combination = dropzeros(sparse(warehouse_combination))
            parcel_evalu = trans * warehouse_combination
            parcel_check = sum(trans; dims = 2)
            parcel = 0
            for i in axes(parcel_evalu, 1)
                for j in axes(parcel_evalu, 2)
                    if parcel_evalu[i, j] == parcel_check[i]
                        parcel += parcels_out[j]
                        break
                    end
                    if j == size(parcel_evalu, 2)
                        if parcel_evalu[i, j] != parcel_check[i]
                            error("Order was not fully dispatched!")
                        end
                    end
                end
            end
        end
        parcel -= size(trans, 1)
    else
        show(X)
        error("Not all SKUs are allocated!")
    end
    return parcel
end

## FLEXIBILITY: average number of warehouses that can individually fulfill each order completely
function FLEXIBILITY(trans::SparseMatrixCSC{Bool,Int64}, X::Array{Bool,2})
    fulfillment = trans * X
    order_sizes = vec(sum(trans; dims = 2))
    flex = 0
    for j in axes(fulfillment, 1)
        for k in axes(fulfillment, 2)
            if fulfillment[j, k] == order_sizes[j]
                flex += 1
            end
        end
    end
    return round(flex / size(trans, 1); digits = 4)
end

# Functions for the random allocation of SKUs to warehouses
## RANDOMALLOCONCE: allocate SKUs randomly in case each SKU can only be assigned once
function RANDOMALLOCONCE(capacity::Vector{<:Real}, sku_weight::Vector{<:Real})
    cap_left = copy(capacity)
    X = zeros(Bool, length(sku_weight), length(capacity))

    # Allocate heaviest SKUs first to prevent blocking, randomize within same weight
    order = sortperm(sku_weight; rev = true)
    for j in order
        valid = findall(k -> cap_left[k] >= sku_weight[j], 1:length(capacity))
        isempty(valid) && error(
            "No warehouse with enough capacity for SKU $j (weight $(sku_weight[j]), max=$(maximum(cap_left)))",
        )
        k = rand(valid)
        X[j, k] = true
        cap_left[k] -= sku_weight[j]
    end

    return X
end

## RANDOMALLOCMULTI: allocate SKUs randomly in case each SKU can be allocated multiple times
function RANDOMALLOCMULTI(
    trans::SparseMatrixCSC{Bool,Int64}, capacity::Vector{<:Real}, sku_weight::Vector{<:Real}
)
    X = RANDOMALLOCONCE(capacity, sku_weight::Vector{<:Real})
    cap_used = vec(sum(X .* sku_weight; dims = 1))
    for d in axes(capacity, 1)
        available = findall(
            i -> !X[i, d] && sku_weight[i] <= capacity[d] - cap_used[d], 1:size(trans, 2)
        )
        while cap_used[d] < capacity[d] && !isempty(available)
            idx = rand(1:length(available))
            i = available[idx]
            if sku_weight[i] <= capacity[d] - cap_used[d]
                X[i, d] = true
                cap_used[d] += sku_weight[i]
                available[idx] = available[end]
                pop!(available)
            else
                available[idx] = available[end]
                pop!(available)
            end
        end
    end
    return X::Array{Bool,2}
end

## RANDOMBENCH: benchmark function to evaluate multiple random allocations
function RANDOMBENCH(
    trans::SparseMatrixCSC{Bool,Int64},
    capacity::Array{Int64,1},
    iterations::Int64,
    sku_weight::Vector{<:Real},
    combination::Array{Array{Array{Int64,1},1},1},
)
    randomcollection = zeros(Int64, iterations)
    flexcollection = zeros(Float64, iterations)
    Threads.@threads for r in 1:iterations
        W = RANDOMALLOCMULTI(trans, capacity, sku_weight)
        randomcollection[r] = PARCELSSEND(trans, W, capacity, combination)
        flexcollection[r] = FLEXIBILITY(trans, W)
    end
    return round(mean(randomcollection)), round(mean(flexcollection); digits = 4)
end

## CHECKCAPACITY: function to test whether the input capacity is viable for the transactional data
function CHECKCAPACITY(capacity::Array{Int64,1}, sku_weight::Vector{<:Real})
    valid = true
    if sum(capacity) < sum(sku_weight)
        error("Abort due to insufficient capacity!")
        valid = false
    end
    return valid::Bool
end
