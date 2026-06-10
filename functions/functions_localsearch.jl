# Minimum number of warehouses needed to cover all SKUs in an order.
# Uses bitmask enumeration (exact) for up to 10 warehouses.
function _min_cover(order_skus, X::Matrix{Bool}, D::Int)
    m = length(order_skus)
    m == 0 && return 0
    full = (UInt32(1) << m) - UInt32(1)
    wh_mask = zeros(UInt32, D)
    @inbounds for (bit_idx, sku) in enumerate(order_skus)
        bit = UInt32(1) << (bit_idx - 1)
        for d in 1:D
            X[sku, d] && (wh_mask[d] |= bit)
        end
    end
    @inbounds for d in 1:D
        wh_mask[d] == full && return 1
    end
    @inbounds for d1 in 1:D-1, d2 in d1+1:D
        (wh_mask[d1] | wh_mask[d2]) == full && return 2
    end
    @inbounds for d1 in 1:D-2, d2 in d1+1:D-1
        m12 = wh_mask[d1] | wh_mask[d2]
        @inbounds for d3 in d2+1:D
            (m12 | wh_mask[d3]) == full && return 3
        end
    end
    # Level 4+: enumerate all 2^D warehouse subsets
    min_p = D
    for ws in UInt16(1):(UInt16(1) << D) - UInt16(1)
        p = count_ones(ws)
        p >= min_p && continue
        cover = UInt32(0)
        for d in 1:D
            if ws & (UInt16(1) << (d - 1)) != 0
                cover |= wh_mask[d]
            end
        end
        cover == full && (min_p = p)
    end
    return min_p
end

# Pair-wise exchange local search on the SKU-warehouse allocation.
# Uses incremental parcel counting and per-warehouse candidate lists.
function LOCALSEARCHCHI!(trans::SparseMatrixCSC{Bool,Int64},
                         X::Matrix{Bool},
                         Q::AbstractMatrix{<:Real},
                         capacity::Vector{Int64},
                         log_results::Bool,
                         ls::Int64,
                         max_ls::Int64,
                         ls_neighborhood::Float64=1.0)
    I = size(X, 1)
    D = size(X, 2)
    J = size(trans, 1)

    # Coappearance ranking (computed once)
    coapp_sum = vec(sum(Q, dims=1))
    coapp_sort = sortperm(coapp_sum, rev=true)

    # Initialize coappearance state matrix
    state = zeros(Float64, I, D)
    CURRENTSTATE!(X, Q, state)

    # Initialize fulfillment tracking: fulfillment[j,d] = # SKUs of order j in warehouse d
    fulfillment = Matrix{Int32}(trans * X)
    replication = sum(capacity) > I
    t_rows = rowvals(trans)

    # Compute initial parcel count
    if !replication
        total_parcels = 0
        @inbounds for j in 1:J, d in 1:D
            fulfillment[j, d] > 0 && (total_parcels += 1)
        end
        total_parcels -= J
        order_parcels = nothing
        trans_t = nothing
    else
        trans_t = SparseMatrixCSC(trans')
        tt_rows = rowvals(trans_t)
        order_parcels = Vector{Int32}(undef, J)
        total_parcels = 0
        @inbounds for j in 1:J
            skus = @view tt_rows[nzrange(trans_t, j)]
            p = _min_cover(skus, X, D)
            order_parcels[j] = p
            total_parcels += p
        end
        total_parcels -= J
    end

    impro_bef = total_parcels

    # Backup structures
    X_backup = similar(X)
    fulfillment_backup = similar(fulfillment)
    order_parcels_backup = replication ? similar(order_parcels) : nothing

    log_results && println("\n  Iter: 0 - parcels: ", impro_bef)

    while ls < max_ls
        ls += 1

        copyto!(X_backup, X)
        copyto!(fulfillment_backup, fulfillment)
        replication && copyto!(order_parcels_backup, order_parcels)
        old_total = total_parcels

        # Search for improvements (returns swapped SKU indices)
        swapped = SEARCHLOOP!(X, Q, coapp_sort, state, ls_neighborhood)

        # Early termination if no swaps were made
        isempty(swapped) && break

        # Incremental parcel update
        affected = replication ? Set{Int}() : nothing
        for sku in unique(swapped)
            @inbounds for d in 1:D
                delta = Int32(X[sku, d]) - Int32(X_backup[sku, d])
                delta == 0 && continue
                @inbounds for idx in nzrange(trans, sku)
                    j = t_rows[idx]
                    old_f = fulfillment[j, d]
                    fulfillment[j, d] += delta
                    if !replication
                        new_f = fulfillment[j, d]
                        old_f > 0 && new_f <= 0 && (total_parcels -= 1)
                        old_f <= 0 && new_f > 0 && (total_parcels += 1)
                    else
                        push!(affected, j)
                    end
                end
            end
        end

        # Replication case: recount min cover for affected orders
        if replication && !isempty(affected)
            tt_rows_ref = rowvals(trans_t)
            for j in affected
                old_p = order_parcels[j]
                skus = @view tt_rows_ref[nzrange(trans_t, j)]
                new_p = _min_cover(skus, X, D)
                total_parcels += (new_p - old_p)
                order_parcels[j] = new_p
            end
        end

        log_results && println("  Iter: ", ls, " - parcels: ", total_parcels)

        if total_parcels < impro_bef
            impro_bef = total_parcels
        else
            copyto!(X, X_backup)
            copyto!(fulfillment, fulfillment_backup)
            replication && copyto!(order_parcels, order_parcels_backup)
            total_parcels = old_total
            break
        end
    end

    return ls
end

function POTENTIAL(state::Matrix{<:Real}, i::Int64, j::Int64, k::Int64, g::Int64)
    state[i, g] - state[i, k] + state[j, k] - state[j, g]
end

function CURRENTSTATE!(X::Matrix{Bool}, Q::AbstractMatrix{<:Real}, state::Matrix{<:Real})
    state .+= Q * X
end

function REFRESHSTATE!(state::Matrix{<:Real},
                       Q::SparseMatrixCSC{<:Real},
                       k::Int64, g::Int64, i::Int64, j::Int64)
    rows = rowvals(Q)
    vals = nonzeros(Q)
    @inbounds for idx in nzrange(Q, j)
        y = rows[idx]
        v = vals[idx]
        state[y, k] += v
        state[y, g] -= v
    end
    @inbounds for idx in nzrange(Q, i)
        y = rows[idx]
        v = vals[idx]
        state[y, k] -= v
        state[y, g] += v
    end
end

# Pair-wise exchange search with per-warehouse candidate lists.
# Returns the list of swapped SKU indices.
function SEARCHLOOP!(X::Matrix{Bool},
                     Q::AbstractMatrix{<:Real},
                     coapp_sort::Vector{Int64},
                     state::Matrix{<:Real},
                     ls_neighborhood::Float64=1.0)
    I = size(X, 1)
    D = size(X, 2)
    K = max(1, ceil(Int, I * ls_neighborhood))
    swapped = Int[]

    # Build per-warehouse candidate lists (ranked by coappearance)
    wh_ranked = [Int[] for _ in 1:D]
    @inbounds for idx in 1:K
        i = coapp_sort[idx]
        for d in 1:D
            X[i, d] && push!(wh_ranked[d], i)
        end
    end

    @fastmath begin
        @inbounds for g in 2:D
            @inbounds for k in 1:D-1
                for i in wh_ranked[k]
                    if X[i, k] && !X[i, g]
                        for j in wh_ranked[g]
                            if X[j, g] && !X[j, k]
                                if POTENTIAL(state, i, j, k, g) > 0
                                    X[i, k] = X[j, g] = false
                                    X[i, g] = X[j, k] = true
                                    REFRESHSTATE!(state, Q, k, g, i, j)
                                    push!(swapped, i)
                                    push!(swapped, j)
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return swapped
end
