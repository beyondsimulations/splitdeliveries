function chi_part(a::Real, b::Real)
    (a-b)^2/b
end

function INDEPENDENT(sum_cond_sku::Vector{<:Real}, J::Int64, i::Int64, j::Int64)
    (Float64(sum_cond_sku[i]) * sum_cond_sku[j])/J
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
                independent = (Float64(sum_cond_sku[row]) * sum_cond_sku[col]) / J
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
                independent = (Float64(sum_cond_sku[i_idx]) * sum_cond_sku[j_idx]) / J
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

## CURVEBALL: margin-preserving randomization of the transaction matrix
## (Strona et al. 2014). Every transaction keeps its size and every SKU its
## frequency, while pairwise dependencies are destroyed. Serves as the
## independence baseline for the dependency gate in CHISQUAREHEUR. Uses a
## private generator so the global random stream is untouched.
function CURVEBALL(trans::SparseMatrixCSC{Bool,Int64}; trades_per_order::Int64 = 10,
                   seed::Int64 = 4242)
    rng = MersenneTwister(seed)
    M, S = size(trans)
    orders = [Int32[] for _ in 1:M]
    rows = rowvals(trans)
    for j in 1:S, idx in nzrange(trans, j)
        push!(orders[rows[idx]], Int32(j))
    end
    # in_a[x] == t (in_b[x] == t): SKU x is in order a (b) of trade t
    in_a = zeros(Int, S)
    in_b = zeros(Int, S)
    # SKUs in only one of the two orders, A's first, and SKUs in both, in A's order
    pool = Int32[]
    shared = Int32[]
    for t in 1:(trades_per_order * M)
        a = rand(rng, 1:M)
        b = rand(rng, 1:M)
        a == b && continue
        A, B = orders[a], orders[b]
        (isempty(A) || isempty(B)) && continue
        for x in A
            in_a[x] = t
        end
        for x in B
            in_b[x] = t
        end
        empty!(pool)
        empty!(shared)
        for x in A
            in_b[x] == t ? push!(shared, x) : push!(pool, x)
        end
        kA = length(pool)
        for x in B
            in_a[x] == t || push!(pool, x)
        end
        (kA == 0 || kA == length(pool)) && continue
        shuffle!(rng, pool)
        # both orders keep their sizes: shared SKUs first, then their part of the pool
        n = length(shared)
        copyto!(A, shared)
        copyto!(A, n + 1, pool, 1, kA)
        copyto!(B, shared)
        copyto!(B, n + 1, pool, kA + 1, length(pool) - kA)
    end
    nnz_total = sum(length, orders)
    I_ = Vector{Int64}(undef, nnz_total)
    J_ = Vector{Int64}(undef, nnz_total)
    k = 0
    for (i, lst) in enumerate(orders), j in lst
        k += 1
        I_[k] = i
        J_[k] = Int64(j)
    end
    return sparse(I_, J_, trues(nnz_total), M, S)
end

## SUPPORTEDDEPMASS: dependency premium carried by flagged SKU pairs that
## coappear in at least two transactions. Pairs observed together only once
## are excluded, since their flagged mass grows with catalog size on
## independent data and carries no diagnostic information.
function SUPPORTEDDEPMASS(dep::SparseMatrixCSC, Q::AbstractMatrix{<:Real})
    rows_d = rowvals(dep)
    vals_d = nonzeros(dep)
    mass = 0.0
    for c in 1:size(dep, 2), idx in nzrange(dep, c)
        r = rows_d[idx]
        r < c || continue
        Q[r, c] >= 2 && (mass += Float64(vals_d[idx]))
    end
    return mass
end

# function to find the bestselling SKUs, borrowed from Catalán and Fisher (2021)
function BESTSELLING_SKUS(capacity_left::Vector{<:Real}, sku_weights::Vector{<:Real})
    (sum(capacity_left)-sum(sku_weights))/(size(capacity_left, 1)-1)
end

# function to allocate the bestselling SKUs to the weights
# sum_nor already holds the per-storage-unit independent coappearances of each
# SKU (the caller divides by the storage size), so it is used as is. The storage
# sizes sku_weights only enter the capacity accounting.
function BESTSELLING_ALLOCATE!(
    sum_nor::Vector{<:Real},
    capacity_left::Vector{<:Real},
    weight::Vector{<:Real},
    sku_weights::Vector{<:Real},
    nor_order::Vector{Int64},
)
    bestselling = BESTSELLING_SKUS(capacity_left, sku_weights)
    normal_weight = Float64.(sum_nor)
    for k in axes(capacity_left, 1)
        if capacity_left[k] > bestselling
            normal_weight .= sum_nor
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
# As in BESTSELLING_ALLOCATE!, sum_nor is already expressed per storage unit and
# sku_weights is only used to book the consumed capacity.
function WHWEIGHT(
    capacity::Vector{Int64}, sum_nor::Vector{<:Real}, sku_weights::Vector{<:Real}
)
    free_capacity::Vector{Float64} = copy(capacity)
    weight = zeros(Float64, size(free_capacity))
    nor_order = sortperm(sum_nor; rev = true)
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
    weight .= weight ./ sum(sum_nor)
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

# Max-heaps per warehouse of the values that FINDDEP compares, entries (value, SKU): dep
# holds the dependent values above zero, nor the independent values above zero. An entry is
# current while its SKU is unallocated and its value equals the value of the current state.
# Every change of a state adds an entry; entries that are not current are dropped at the
# top. nan is set once a value is NaN; FINDDEP decides from then on.
mutable struct DEPQUEUE{S<:Real,W<:Real}
    state_dep::Matrix{S}
    state_nor::Matrix{S}
    sku_weight::Vector{W}
    allocated::Vector{Bool}
    dep::Vector{Vector{Tuple{Float64,Int64}}}
    nor::Vector{Vector{Tuple{Float64,Int64}}}
    aside::Vector{Tuple{Float64,Int64}}
    nan::Bool
end
function DEPQUEUE(
    state_dep::Matrix{<:Real},
    state_nor::Matrix{<:Real},
    sku_weight::Vector{<:Real},
    allocated::Vector{Bool},
)
    K = size(state_dep, 2)
    return DEPQUEUE(
        state_dep,
        state_nor,
        sku_weight,
        allocated,
        [Tuple{Float64,Int64}[] for _ in 1:K],
        [Tuple{Float64,Int64}[] for _ in 1:K],
        Tuple{Float64,Int64}[],
        false,
    )
end

# dependent and independent value of SKU j for warehouse k, computed as in FINDDEP
@inline function DEPVALUES(queue::DEPQUEUE, j::Integer, k::Int64)
    weight_inv = 1.0 / queue.sku_weight[j]
    pot_dep = queue.state_dep[j, k] * weight_inv
    pot_nor = queue.state_nor[j, k] * weight_inv
    if pot_dep > 0
        pot_dep += pot_nor
    end
    return pot_dep, pot_nor
end

# entry e of the heap of warehouse k with the dependent (v = 1) or the independent (v = 2)
# values is current
CURRENT(queue::DEPQUEUE, e::Tuple{Float64,Int64}, v::Int64, k::Int64) =
    !queue.allocated[e[2]] && DEPVALUES(queue, e[2], k)[v] == e[1]

# heap order: higher value first, ties to the smaller SKU
HEAPBEFORE(a::Tuple{Float64,Int64}, b::Tuple{Float64,Int64}) =
    a[1] > b[1] || (a[1] == b[1] && a[2] < b[2])

function HEAPPUSH!(h::Vector{Tuple{Float64,Int64}}, e::Tuple{Float64,Int64})
    push!(h, e)
    c = length(h)
    @inbounds while c > 1 && HEAPBEFORE(h[c], h[c >> 1])
        h[c], h[c >> 1] = h[c >> 1], h[c]
        c >>= 1
    end
end

function HEAPPOP!(h::Vector{Tuple{Float64,Int64}})
    top = h[1]
    e = pop!(h)
    n = length(h)
    if n > 0
        p = 1
        @inbounds while 2p <= n
            c = 2p < n && HEAPBEFORE(h[2p + 1], h[2p]) ? 2p + 1 : 2p
            HEAPBEFORE(h[c], e) || break
            h[p] = h[c]
            p = c
        end
        h[p] = e
    end
    return top
end

# add the values of the unallocated SKU j for warehouse k after a change of its state; a
# heap with more than twice as many entries as SKUs keeps only its current entries, each
# once
function PUSHDEP!(queue::DEPQUEUE, j::Integer, k::Int64)
    pot_dep, pot_nor = DEPVALUES(queue, j, k)
    queue.nan |= isnan(pot_dep) | isnan(pot_nor)
    pot_dep > 0 && HEAPPUSH!(queue.dep[k], (pot_dep, Int64(j)))
    pot_nor > 0 && HEAPPUSH!(queue.nor[k], (pot_nor, Int64(j)))
    for (h, v) in ((queue.dep[k], 1), (queue.nor[k], 2))
        if length(h) > 2 * length(queue.allocated)
            filter!(e -> CURRENT(queue, e, v, k), h)
            unique!(sort!(h; lt = HEAPBEFORE))
        end
    end
end

# first current entry of the heap of warehouse k with the dependent (v = 1) or independent
# (v = 2) values whose SKU fits into k, or (0.0, 0). Entries that are not current are
# dropped, current entries of SKUs that do not fit stay in the heap.
function TOPFIT!(queue::DEPQUEUE, v::Int64, k::Int64, capacity_left::Vector{<:Real})
    h = v == 1 ? queue.dep[k] : queue.nor[k]
    top = (0.0, 0)
    while !isempty(h)
        e = h[1]
        if !CURRENT(queue, e, v, k)
            HEAPPOP!(h)
        elseif capacity_left[k] < queue.sku_weight[e[2]]
            push!(queue.aside, HEAPPOP!(h))
        else
            top = e
            break
        end
    end
    for e in queue.aside
        HEAPPUSH!(h, e)
    end
    empty!(queue.aside)
    return top
end

# FINDDEP from the heaps of warehouse k: the highest dependent value among the unallocated
# SKUs that fit into k, the first SKU with this value and the highest independent value
# among these SKUs. A value of zero stands for none above zero; ADDDEPENDENT! decides as
# with FINDDEP.
function QUEUEDEP(queue::DEPQUEUE, k::Int64, capacity_left::Vector{<:Real})
    max_dep_val, i = TOPFIT!(queue, 1, k, capacity_left)
    i == 0 && return 0.0, 1, 0.0
    return max_dep_val, i, TOPFIT!(queue, 2, k, capacity_left)[1]
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
    queue::Union{DEPQUEUE,Nothing},
)
    while capacity_left[k] > 0 && n_allocated[] < size(X, 1)
        max_dep_val, i, max_nor_val = if queue === nothing || queue.nan
            FINDDEP(X, k, state_dep, state_nor, allocated, sku_weight, capacity_left)
        else
            QUEUEDEP(queue, k, capacity_left)
        end
        max_dep_val <= 0 && break
        if max_dep_val >= max_nor_val && capacity_left[k] >= sku_weight[i]
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
                queue,
            )
        else
            break
        end
    end
end

# function to check the dependencies to already allocated SKUs: returns the highest
# dependent value among the unallocated SKUs that fit into warehouse k, the first SKU
# with this value and the highest independent value. The dependent value is NaN if any
# value is NaN.
function FINDDEP(
    X::Array{Bool,2},
    k::Int64,
    state_dep::Matrix{<:Real},
    state_nor::Matrix{<:Real},
    allocated::Vector{Bool},
    sku_weight::Vector{<:Real},
    capacity_left::Vector{<:Real},
)
    max_dep_val, i, max_nor_val = -Inf, 1, -Inf
    nan = false
    @inbounds for j in 1:size(X, 1)
        pot_dep = pot_nor = 0.0
        if allocated[j] == 0 && capacity_left[k] >= sku_weight[j]
            weight_inv = 1.0 / sku_weight[j]
            pot_dep = state_dep[j, k] * weight_inv
            pot_nor = state_nor[j, k] * weight_inv
            if pot_dep > 0
                pot_dep += pot_nor
            end
        end
        nan |= isnan(pot_dep) | isnan(pot_nor)
        if pot_dep > max_dep_val
            max_dep_val, i = pot_dep, j
        end
        if pot_nor > max_nor_val
            max_nor_val = pot_nor
        end
    end
    nan && (max_dep_val = NaN)
    return max_dep_val, i, max_nor_val
end

# function to allocate a selcted product to a selected warehouse; every changed state of an
# unallocated SKU is added to the heaps of queue
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
    queue::Union{DEPQUEUE,Nothing} = nothing,
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
                            if queue !== nothing
                                PUSHDEP!(queue, m, d_old[ei])
                                PUSHDEP!(queue, m, d_new[ei])
                            end
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
                        queue === nothing || PUSHDEP!(queue, m, d_placed)
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
                queue === nothing || PUSHDEP!(queue, j, k)
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
    n_unallocated = count(unallocated)

    best_allocation = zeros(Float64, size(Q, 1))
    q_rows = rowvals(Q)
    q_vals = nonzeros(Q)

    # With positive finite storage requirements and finite coappearances of at least zero,
    # the states of a warehouse only grow while it is filled and the next copy comes from
    # a max-heap of these states (FILLTOP!); otherwise each copy scans all SKUs
    use_heap =
        all(w -> 0 < w < Inf && 1.0 / w < Inf, sku_weight) &&
        all(v -> isfinite(v) && v >= 0, q_vals)
    heap = Tuple{Float64,Int64}[]

    @inbounds for d in 1:size(capacity_left, 1)
        # heap of warehouse d: (state, SKU) of the SKUs not in d with a state above zero
        if use_heap && capacity_left[d] >= min_weight
            empty!(heap)
            for i in 1:size(X, 1)
                if !X[i, d] && state[i, d] > 0
                    push!(heap, (state[i, d], i))
                end
            end
            sort!(heap; lt = HEAPBEFORE)
        end
        while capacity_left[d] >= min_weight
            if !use_heap
                fill!(best_allocation, 0.0)

                @inbounds @simd for i in 1:size(Q, 1)
                    if !X[i, d] && capacity_left[d] >= sku_weight[i]
                        best_allocation[i] = state[i, d]
                    end
                end
            end

            # Prioritize SKUs that have no allocation yet
            fallback = 0
            if n_unallocated > 0
                for i in 1:size(X, 1)
                    if unallocated[i] && capacity_left[d] >= sku_weight[i]
                        fallback = i
                        break
                    end
                end
            end
            if fallback > 0
                X[fallback, d] = 1
                capacity_left[d] -= sku_weight[fallback]
                unallocated[fallback] = false
                n_unallocated -= 1
            else
                best = if use_heap
                    FILLTOP!(heap, X, state, capacity_left, sku_weight, d)
                else
                    findmax(best_allocation)
                end
                if best[1] > 0
                    X[best[2], d] = 1
                    capacity_left[d] -= sku_weight[best[2]]

                    @inbounds for idx in nzrange(Q, best[2])
                        i_row = q_rows[idx]
                        if sku_weight[i_row] > 0
                            value = state[i_row, d]
                            state[i_row, d] += q_vals[idx] / sku_weight[i_row]
                            if use_heap && !X[i_row, d] && state[i_row, d] != value
                                HEAPPUSH!(heap, (state[i_row, d], Int64(i_row)))
                            end
                        end
                    end
                else
                    capacity_left[d] = 0
                    break
                end
            end
        end
    end
end

# first entry of the heap h of warehouse d whose SKU is not in d, whose value equals the
# state of the SKU and whose SKU fits into d, or (0.0, 0) if there is none; the entries
# before it are dropped. While d is filled, X[:, d] is only set, the states of d only grow
# and the capacity left of d only falls, so a dropped entry never qualifies again.
function FILLTOP!(
    h::Vector{Tuple{Float64,Int64}},
    X::Array{Bool,2},
    state::Matrix{Float64},
    capacity_left::Vector{<:Real},
    sku_weight::Vector{<:Real},
    d::Int64,
)
    while !isempty(h)
        e = HEAPPOP!(h)
        if !X[e[2], d] && state[e[2], d] == e[1] && capacity_left[d] >= sku_weight[e[2]]
            return e
        end
    end
    return (0.0, 0)
end

# Capacities for the allocation of the single copies: the buffer sum(capacity) - target is
# held free in equal whole shares in every warehouse but the largest (capacity sorted
# decreasingly). A share is capped at the warehouse's capacity; what cannot be held free
# there is held free in the largest warehouse.
function RESERVE(capacity::Vector{Int64}, target::Int64)
    K = length(capacity)
    R = sum(capacity) - target
    r = zeros(Int64, K)
    while R > 0
        open = [k for k in 2:K if r[k] < capacity[k]]
        if isempty(open)
            r[1] += R
            break
        end
        q = max(R ÷ length(open), 1)
        for k in open
            a = min(q, capacity[k] - r[k], R)
            r[k] += a
            R -= a
        end
    end
    return capacity .- r
end

# Fill the capacity left in each warehouse, one warehouse after the other. A copy of SKU i
# in warehouse d is valued by the coappearances of i with the SKUs in d that share no
# warehouse with i, per storage unit; ties go to the higher order frequency. If use_pairs is
# false or no candidate has a positive value, the next SKU of nor_order that is not in d and
# fits is copied.
function FILLRESERVE!(
    X::Matrix{Bool},
    Q::SparseMatrixCSC{<:Real},
    cap_left::Vector{<:Real},
    sku_weight::Vector{<:Real},
    ordered_skus::Vector{<:Real},
    nor_order::Vector{Int64},
    use_pairs::Bool,
)
    I, D = size(X)
    # mask[i]: one bit per warehouse that holds SKU i
    mask = [reduce(|, (UInt64(X[i, d]) << (d - 1) for d in 1:D)) for i in 1:I]
    # P[i, d]: coappearances of i with the SKUs in d that share no warehouse with i
    P = zeros(Float64, I, D)
    rows = rowvals(Q)
    vals = nonzeros(Q)
    if use_pairs
        for i in 1:I, idx in nzrange(Q, i)
            j = rows[idx]
            (j == i || mask[i] & mask[j] != 0) && continue
            for d in 1:D
                X[j, d] && (P[i, d] += vals[idx])
            end
        end
    end
    minw = minimum(sku_weight)
    for d in 1:D
        while cap_left[d] >= minw
            best, bf, bi = 0.0, -1.0, 0
            if use_pairs
                @inbounds for i in 1:I
                    (X[i, d] || cap_left[d] < sku_weight[i]) && continue
                    v = P[i, d] / sku_weight[i]
                    if v > best || (v == best && v > 0 && ordered_skus[i] > bf)
                        best, bf, bi = v, Float64(ordered_skus[i]), i
                    end
                end
            end
            if bi == 0
                for i in nor_order
                    if !X[i, d] && cap_left[d] >= sku_weight[i]
                        bi = i
                        break
                    end
                end
            end
            bi == 0 && break
            if use_pairs
                # pairs of bi that become stored together stop counting in every warehouse;
                # SKUs that are not in d can be brought together with bi there
                for idx in nzrange(Q, bi)
                    j = rows[idx]
                    (j == bi || mask[bi] & mask[j] != 0) && continue
                    if X[j, d]
                        for e in 1:D
                            e == d && continue
                            X[j, e] && (P[bi, e] -= vals[idx])
                            X[bi, e] && (P[j, e] -= vals[idx])
                        end
                    else
                        P[j, d] += vals[idx]
                    end
                end
            end
            X[bi, d] = true
            mask[bi] |= UInt64(1) << (d - 1)
            cap_left[d] -= sku_weight[bi]
        end
    end
    return nothing
end
