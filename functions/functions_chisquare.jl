function chi_part(a::Real, b::Real)
    (a-b)^2/b
end

function INDEPENDENT(sum_cond_sku::Vector{<:Real}, J::Int64, i::Int64, j::Int64)
    (sum_cond_sku[i] * sum_cond_sku[j])/J
end

# function used to perform the chi square test of independence upon
# the coappearance matrix Q
function HYOPTHESISTEST_SPARSE(
    Q::AbstractMatrix{<:Real},
    I::Int64,
    J::Int64,
    sig::Float64,
    sum_cond_sku::Vector{<:Real},
    min_effect::Float64 = 0.0,
)
    M = convert(Int64, ((I^2) - I) / 2)
    accept = cquantile(Chisq(1), sig / M)

    if Q isa SparseMatrixCSC
        q_rows = rowvals(Q)
        q_vals = nonzeros(Q)
        # Per-column output arrays (no sharing between threads)
        col_rows = [Int32[] for _ in 1:I]
        col_cols = [Int32[] for _ in 1:I]
        col_vals = [Float32[] for _ in 1:I]

        Threads.@threads for col in 1:I
            @inbounds for idx in nzrange(Q, col)
                row = q_rows[idx]
                row >= col && continue
                qij = q_vals[idx]
                independent = (sum_cond_sku[row] * sum_cond_sku[col]) / J
                qij <= independent && continue

                chi_nr = J - sum_cond_sku[row]
                chi_nd = J - sum_cond_sku[col]
                chi_yn = sum_cond_sku[row] - qij
                chi_ny = sum_cond_sku[col] - qij
                chi_nn = chi_nd - chi_yn

                ind_yn = (sum_cond_sku[row] * chi_nd) / J
                ind_ny = (sum_cond_sku[col] * chi_nr) / J
                ind_nn = (chi_nr * chi_nd) / J

                chi =
                    chi_part(qij, independent) +
                    chi_part(chi_yn, ind_yn) +
                    chi_part(chi_ny, ind_ny) +
                    chi_part(chi_nn, ind_nn)

                if chi > accept
                    val = qij - independent
                    if val >= independent * min_effect
                        push!(col_rows[col], row);
                        push!(col_cols[col], col);
                        push!(col_vals[col], val)
                        push!(col_rows[col], col);
                        push!(col_cols[col], row);
                        push!(col_vals[col], val)
                    end
                end
            end
        end

        # Merge per-column results
        dep_rows = reduce(vcat, col_rows)
        dep_cols = reduce(vcat, col_cols)
        dep_vals = reduce(vcat, col_vals)
    else
        dep_rows = Int32[]
        dep_cols = Int32[]
        dep_vals = Float32[]
        @inbounds for j_idx in 2:I
            @inbounds for i_idx in 1:(j_idx - 1)
                independent = (sum_cond_sku[i_idx] * sum_cond_sku[j_idx]) / J
                qij = Q[i_idx, j_idx]
                qij <= independent && continue

                chi_nr = J - sum_cond_sku[i_idx]
                chi_nd = J - sum_cond_sku[j_idx]
                chi_yn = sum_cond_sku[i_idx] - qij
                chi_ny = sum_cond_sku[j_idx] - qij
                chi_nn = chi_nd - chi_yn

                ind_yn = (sum_cond_sku[i_idx] * chi_nd) / J
                ind_ny = (sum_cond_sku[j_idx] * chi_nr) / J
                ind_nn = (chi_nr * chi_nd) / J

                chi =
                    chi_part(qij, independent) +
                    chi_part(chi_yn, ind_yn) +
                    chi_part(chi_ny, ind_ny) +
                    chi_part(chi_nn, ind_nn)

                if chi > accept
                    val = qij - independent
                    if val >= independent * min_effect
                        push!(dep_rows, i_idx);
                        push!(dep_cols, j_idx);
                        push!(dep_vals, val)
                        push!(dep_rows, j_idx);
                        push!(dep_cols, i_idx);
                        push!(dep_vals, val)
                    end
                end
            end
        end
    end

    return sparse(dep_rows, dep_cols, dep_vals, I, I)
end

# function to find the bestselling SKUs, borrowed from Catalán and Fisher (2021)
function BESTSELLING_SKUS(capacity_left::Vector{<:Real}, sku_weights::Vector{<:Real})
    (sum(capacity_left)-sum(sku_weights))/(size(capacity_left, 1)-1)
end

# function to pre-replicate the bestselling SKUs to all warehouses
# before the Phase 1 dependency-aware allocation loop.
# Weight-aware: budget and capacity checks use sku_weight per SKU.
function REPLICATEALL!(
    X::Matrix{Bool},
    dep::SparseMatrixCSC{<:Real},
    Q::SparseMatrixCSC{<:Real},
    sum_dep::Vector{<:Real},
    sum_nor::Vector{<:Real},
    state_dep::Matrix{<:Real},
    state_nor::Matrix{<:Real},
    cap_left::Vector{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    nor_order::Vector{Int64},
    n_allocated::Ref{Int},
)
    D = size(X, 2)
    budget = BESTSELLING_SKUS(cap_left, sku_weight)
    if budget <= 0 || D <= 1
        return nothing
    end

    budget_left = budget
    q_rows = rowvals(Q)
    q_vals = nonzeros(Q)
    nonuniform = !all(w -> w == sku_weight[1], sku_weight)

    # Pre-build remaining weights vector for feasibility checks (non-uniform only)
    if nonuniform
        remaining_w = [sku_weight[j] for j in 1:size(X, 1) if !allocated[j]]
        sort!(remaining_w; rev = true)
    end

    for pos in 1:length(nor_order)
        i = nor_order[pos]
        w_i = sku_weight[i]

        # Stop if this SKU exceeds remaining budget
        budget_left - w_i < 0 && break

        # Check that ALL warehouses have enough capacity
        can_fit = true
        @inbounds for k in 1:D
            if cap_left[k] < w_i
                can_fit = false
                break
            end
        end
        can_fit || continue

        # Check feasibility before committing (non-uniform weights only)
        if nonuniform
            @inbounds for k in 1:D
                cap_left[k] -= w_i
            end
            # Remove this SKU's weight from remaining
            idx_rm = searchsortedfirst(remaining_w, w_i; rev = true)
            deleteat!(remaining_w, idx_rm)
            if !FEASIBLE_REMAINING(cap_left, remaining_w)
                @inbounds for k in 1:D
                    cap_left[k] += w_i
                end
                # Re-insert the weight
                insert!(remaining_w, idx_rm, w_i)
                break
            end
        else
            @inbounds for k in 1:D
                cap_left[k] -= w_i
            end
        end

        # Commit allocation
        @inbounds for k in 1:D
            X[i, k] = true
            state_dep[i, k] = 0.0
            state_nor[i, k] = 0.0
        end

        allocated[i] = true
        n_allocated[] += 1
        sum_dep[i] = 0.0
        sum_nor[i] = 0.0

        # Update state for all unallocated SKUs in all warehouses
        @inbounds for idx in nzrange(Q, i)
            j = q_rows[idx]
            if !allocated[j]
                dv = dep[j, i]
                nv = q_vals[idx] - dv
                for k in 1:D
                    state_dep[j, k] += dv
                    state_nor[j, k] += nv
                end
            end
        end

        budget_left -= w_i
    end
end

# function to allocate the bestselling SKUs to the weights
function BESTSELLING_ALLOCATE!(
    sum_nor::Vector{<:Real},
    capacity_left::Vector{<:Real},
    weight::Vector{<:Real},
    sku_weights::Vector{<:Real},
    nor_order::Vector{Int64},
)
    bestselling = BESTSELLING_SKUS(capacity_left, sku_weights)
    normal_weight = sum_nor ./ sku_weights
    for k in axes(capacity_left, 1)
        if capacity_left[k] > bestselling
            normal_weight .= sum_nor ./ sku_weights
            bestselling_space_left = copy(bestselling)
            pos = 1
            while pos <= length(nor_order) &&
                bestselling_space_left - sku_weights[nor_order[pos]] >= 0
                i = nor_order[pos]
                weight[k] += normal_weight[i]
                capacity_left[k] -= sku_weights[i]
                bestselling_space_left -= sku_weights[i]
                normal_weight[i] = 0.0
                pos += 1
            end
        end
    end
    return normal_weight
end

# function to calculate the weight of each warehouse. It shows us the density of 
# the independent coappearances in each warehouse if we were to allocate
# all SKUs according to the highest independent coappearances.
function WHWEIGHT(
    capacity::Vector{Int64}, sum_nor::Vector{<:Real}, sku_weights::Vector{<:Real}
)
    free_capacity::Vector{Float64} = copy(capacity)
    weight = zeros(Float64, size(free_capacity))
    nor_order = sortperm(sum_nor ./ sku_weights; rev = true)
    normal = BESTSELLING_ALLOCATE!(sum_nor, free_capacity, weight, sku_weights, nor_order)
    normal_weight = copy(normal)
    pos = 1
    for k in 1:size(capacity, 1)
        while pos <= length(nor_order)
            i = nor_order[pos]
            if normal_weight[i] <= 0
                pos += 1
                continue
            end
            if free_capacity[k] < sku_weights[i]
                break
            end
            weight[k] += normal_weight[i]
            free_capacity[k] -= sku_weights[i]
            normal_weight[i] = 0
            pos += 1
        end
    end
    weight .= weight ./ sum(sum_nor ./ sku_weights)
    return weight::Vector{Float64}
end

# function to find the largest warehouse that still has
# space left
function WHSPACE(capacity_left::Vector{<:Real}, space::Int64)
    k_ind = 0
    for k in 1:size(capacity_left, 1)
        if capacity_left[k] >= space
            k_ind = k
            break
        end
    end
    if k_ind == 0
        k_ind = argmax(capacity_left)
    end
    return k_ind::Int64
end

# function to select the best SKU for an allocation
function SELECTIK(
    sum_dep::Vector{<:Real},
    sum_nor::Vector{<:Real},
    weight::Vector{Float64},
    capacity_left::Vector{<:Real},
    state_dep::Matrix{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    nor_order::Vector{Int64},
    nor_pos::Ref{Int},
)
    while nor_pos[] <= length(nor_order) && allocated[nor_order[nor_pos[]]]
        nor_pos[] += 1
    end
    if nor_pos[] <= length(nor_order)
        i = nor_order[nor_pos[]]
    else
        i = 1
        for j in 1:length(allocated)
            if !allocated[j]
                i = j
                break
            end
        end
    end
    k_ind = WHSPACE(capacity_left, ceil(Int64, sku_weight[i]))
    k_max = findmax(capacity_left)[2]
    pot_dep = WHPOTDEP(capacity_left, i, state_dep, sku_weight)
    k_dep = findmax(pot_dep)[2]
    if pot_dep[k_dep] > 0 &&
        k_ind != k_dep &&
        sum_dep[i] + sum_nor[i] * weight[k_dep] > sum_nor[i] * weight[k_ind]
        k = k_dep
    elseif k_ind != k_max &&
        sum_dep[i] + sum_nor[i] * weight[k_max] > sum_nor[i] * weight[k_ind]
        k = k_max
    else
        k = k_ind
    end
    return i::Int64, k::Int64
end

# function to check for each warehouse with free space whether
# SKU i has significant dependencies to other already allocated SKUs
function WHPOTDEP(
    capacity_left::Vector{<:Real},
    i::Int64,
    state_dep::Matrix{<:Real},
    sku_weight::Vector{<:Real},
)
    pot_dep = zeros(Float64, size(capacity_left, 1))
    w_inv = 1.0 / sku_weight[i]
    @inbounds for k in 1:size(capacity_left, 1)
        if capacity_left[k] >= sku_weight[i]
            pot_dep[k] = state_dep[i, k] * w_inv
        end
    end
    return pot_dep::Vector{Float64}
end

# function to remove every so far assigned dependent SKU-pair from
# the coappearance matrix dep to prevent the allocation bias described
# in our article.
function REMOVEALLOC(
    X::Array{Bool,2}, Q::SparseMatrixCSC{<:Real}, dep::SparseMatrixCSC{<:Real}
)
    q_rows = rowvals(Q)
    q_vals = nonzeros(Q)
    d_rows = rowvals(dep)
    d_vals = nonzeros(dep)

    # Allocate only new values array, reuse Q's sparsity structure
    new_vals = Vector{Float32}(undef, length(q_vals))

    for col in 1:size(Q, 2)
        d_range = nzrange(dep, col)
        d_pos = first(d_range)
        d_end = last(d_range)

        for q_idx in nzrange(Q, col)
            row = q_rows[q_idx]
            qval = Float64(q_vals[q_idx])

            # Merge-scan dep for matching row
            while d_pos <= d_end && d_rows[d_pos] < row
                d_pos += 1
            end
            dval = (d_pos <= d_end && d_rows[d_pos] == row) ? Float64(d_vals[d_pos]) : 0.0

            if dval > 0
                # a pair stored together in ANY warehouse has realised its
                # dependency premium; credit only the baseline coappearance
                colocated = false
                for k in 1:size(X, 2)
                    if X[row, k] && X[col, k]
                        colocated = true
                        break
                    end
                end
                new_vals[q_idx] = Float32(colocated ? qval - dval : qval)
            else
                new_vals[q_idx] = Float32(qval)
            end
        end
    end

    return SparseMatrixCSC(size(Q, 1), size(Q, 2), Q.colptr, q_rows, new_vals)
end

# function to check for all unallocated SKUs whether they have positive 
# dependencies to the SKUs in the warehouse k the last SKU was allocated 
# to. If so, check whether the dependencies are expected to dominate the 
# independent coapperances. If yes, allocate the corresponding SKUs to 
# the warehouse k.
function ADDDEPENDENT!(
    X::Matrix{Bool},
    Q::AbstractMatrix{<:Real},
    capacity_left::Vector{<:Real},
    k::Int64,
    dep::AbstractMatrix{<:Real},
    sum_dep::Vector{<:Real},
    sum_nor::Vector{<:Real},
    state_dep::Matrix{<:Real},
    state_nor::Matrix{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    n_allocated::Ref{Int},
)
    pot_dep = zeros(Float64, size(X, 1))
    pot_nor = zeros(Float64, size(X, 1))
    while capacity_left[k] > 0 && n_allocated[] < size(X, 1)
        FINDDEP!(
            X,
            k,
            state_dep,
            state_nor,
            pot_dep,
            pot_nor,
            allocated,
            sku_weight,
            capacity_left,
        )
        max_dep_val, i = findmax(pot_dep)
        max_dep_val <= 0 && break
        max_nor_val = findmax(pot_nor)[1]
        if pot_dep[i] >= max_nor_val && capacity_left[k] >= sku_weight[i]
            ALLOCATEONE!(
                X,
                dep,
                Q,
                sum_dep,
                sum_nor,
                state_dep,
                state_nor,
                capacity_left,
                allocated,
                sku_weight,
                i,
                k,
                n_allocated,
            )
        else
            break
        end
    end
end

# function to check the dependencies to already allocated SKUs
function FINDDEP!(
    X::Array{Bool,2},
    k::Int64,
    state_dep::Matrix{<:Real},
    state_nor::Matrix{<:Real},
    pot_dep::Vector{<:Real},
    pot_nor::Vector{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    capacity_left::Vector{<:Real},
)
    @inbounds @simd for j in 1:size(X, 1)
        if allocated[j] == 0 && capacity_left[k] >= sku_weight[j]
            weight_inv = 1.0 / sku_weight[j]
            pot_dep[j] = state_dep[j, k] * weight_inv
            pot_nor[j] = state_nor[j, k] * weight_inv
            if pot_dep[j] > 0
                pot_dep[j] += pot_nor[j]
            end
        else
            pot_dep[j] = pot_nor[j] = 0
        end
    end
end

# function to allocate a selcted product to a selected warehouse
function ALLOCATEONE!(
    X::Array{Bool,2},
    dep::AbstractMatrix{<:Real},
    Q::SparseMatrixCSC{<:Real},
    sum_dep::Vector{<:Real},
    sum_nor::Vector{<:Real},
    state_dep::Matrix{<:Real},
    state_nor::Matrix{<:Real},
    capacity_left::Vector{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    i::Int64,
    k::Int64,
    n_allocated::Ref{Int},
)
    if !allocated[i]
        if capacity_left[k] < sku_weight[i]
            k_safe = argmax(capacity_left)
            if capacity_left[k_safe] >= sku_weight[i]
                k = k_safe
            elseif !all(w -> w == sku_weight[1], sku_weight)
                d_placed, evicted, d_old, d_new = SWAPREPAIR!(
                    X, capacity_left, sku_weight, i
                )
                # Repair state arrays for each moved SKU (sparse Q)
                q_rows = rowvals(Q)
                q_vals = nonzeros(Q)
                for (ei, j) in enumerate(evicted)
                    @inbounds for idx in nzrange(Q, j)
                        m = q_rows[idx]
                        if !allocated[m]
                            dv = dep[m, j]
                            nv = q_vals[idx] - dv
                            state_dep[m, d_old[ei]] -= dv
                            state_nor[m, d_old[ei]] -= nv
                            state_dep[m, d_new[ei]] += dv
                            state_nor[m, d_new[ei]] += nv
                        end
                    end
                end
                # credit SKU i's own coappearances to its new warehouse so
                # later dependency chaining sees the placement
                @inbounds for idx in nzrange(Q, i)
                    m = q_rows[idx]
                    if !allocated[m]
                        dv = dep[m, i]
                        state_dep[m, d_placed] += dv
                        state_nor[m, d_placed] += q_vals[idx] - dv
                    end
                end
                allocated[i] = true
                n_allocated[] += 1
                sum_dep[i] = 0
                sum_nor[i] = 0
                return nothing
            else
                error(
                    "Cannot allocate SKU $i (weight $(sku_weight[i])): no warehouse has enough capacity.",
                )
            end
        end
        X[i, k] = 1
        sum_dep[i] = 0
        sum_nor[i] = 0
        capacity_left[k] -= sku_weight[i]
        allocated[i] = true
        n_allocated[] += 1
        state_dep[i, k] = 0
        state_nor[i, k] = 0
        q_rows = rowvals(Q)
        q_vals = nonzeros(Q)
        @inbounds for idx in nzrange(Q, i)
            j = q_rows[idx]
            if !allocated[j]
                dv = dep[j, i]
                state_dep[j, k] += dv
                state_nor[j, k] += q_vals[idx] - dv
            end
        end
    else
        error("This product is already allocated!")
    end
end

## check whether all warehouse except the last one are already full. If
## that is the case just allocate the remaining SKUs yet not allocated there.
function FILLLAST!(
    X::Array{Bool,2},
    capacity_left::Vector{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    n_allocated::Ref{Int},
)
    if sum(capacity_left[1:(size(capacity_left, 1) - 1)]) <= 0
        last_d = size(capacity_left, 1)
        @inbounds for i in 1:length(allocated)
            if allocated[i] == false
                if capacity_left[last_d] < sku_weight[i]
                    k_safe = argmax(capacity_left)
                    if capacity_left[k_safe] >= sku_weight[i]
                        last_d = k_safe
                    else
                        SWAPREPAIR!(X, capacity_left, sku_weight, i)
                        allocated[i] = true
                        n_allocated[] += 1
                        continue
                    end
                end
                X[i, last_d] = 1
                capacity_left[last_d] -= sku_weight[i]
                allocated[i] = true
                n_allocated[] += 1
            end
        end
    end
end

# function to allocate the SKUs with the highest potential allocation 
# value to each warehouse with leftover storage space until it is full
function FILLUP!(
    X::Array{Bool,2},
    Q::AbstractMatrix{<:Real},
    capacity_left::Vector{<:Real},
    sku_weight::Vector{<:Real},
)
    min_weight = minimum(sku_weight)

    state = Matrix{Float64}(Q * X)
    @inbounds for i in 1:size(X, 1)
        if sku_weight[i] > 0
            w_inv = 1.0 / sku_weight[i]
            @inbounds @simd for d in 1:size(X, 2)
                state[i, d] *= w_inv
            end
        end
    end

    # Track which SKUs still need their first allocation
    unallocated = [sum(@view(X[i, :])) == 0 for i in 1:size(X, 1)]

    best_allocation = zeros(Float64, size(Q, 1))
    q_rows = rowvals(Q)
    q_vals = nonzeros(Q)

    @inbounds for d in 1:size(capacity_left, 1)
        while capacity_left[d] >= min_weight
            fill!(best_allocation, 0.0)

            @inbounds @simd for i in 1:size(Q, 1)
                if !X[i, d] && capacity_left[d] >= sku_weight[i]
                    best_allocation[i] = state[i, d]
                end
            end

            # Prioritize SKUs that have no allocation yet
            fallback = 0
            for i in 1:size(X, 1)
                if unallocated[i] && capacity_left[d] >= sku_weight[i]
                    fallback = i
                    break
                end
            end
            if fallback > 0
                X[fallback, d] = 1
                capacity_left[d] -= sku_weight[fallback]
                unallocated[fallback] = false
            else
                best = findmax(best_allocation)
                if best[1] > 0
                    X[best[2], d] = 1
                    capacity_left[d] -= sku_weight[best[2]]

                    weight_inv = sku_weight[best[2]] > 0 ? 1.0 / sku_weight[best[2]] : 0.0
                    @inbounds for idx in nzrange(Q, best[2])
                        state[q_rows[idx], d] += q_vals[idx] * weight_inv
                    end
                else
                    capacity_left[d] = 0
                    break
                end
            end
        end
    end
end
