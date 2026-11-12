function CHISQUAREHEUR(
    trans::SparseMatrixCSC{Bool,Int64},
    capacity::Vector{Int64},
    sig::Float64,
    max_ls::Int64,
    sku_weight::Vector{<:Real},
    log_results::Bool,
    min_effect::Float64 = 0.0,
    ls_neighborhood::Float64 = 1.0,
    gate_ratio::Float64 = 2.0,
)
    # Number of transactions
    J = size(trans, 1)
    # Number of products
    I = size(trans, 2)
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity; rev = true)
    # Set if the independence gate removes the dependency premiums
    gate_fired = false
    # Create Coapperance Matrix
    log_results == true ? print("\n  creating coappearance matrix.") : nothing
    Q = COAPPEARENCE(trans, sku_weight)
    # Start the main part
    if CHECKCAPACITY(capacity, sku_weight)
        # Create Vector with number of transactions containing each SKU
        ordered_skus = [Q[i, i] for i in axes(Q, 1)]
        # Clean the principle diagonal
        CLEANPRINCIPLE!(Q)
        # Chi Square Test
        ## We start the chi-square heuristic by performing chi-square tests
        ## for all unique SKU-combinations. For each undirected SKU-pair we
        ## evaluate a chi-square test and either accept or reject the null
        ## hypothesis of independence. Our goal is a matrix dep with the
        ## “net-benefit” of all unique positive dependent SKU-combinations
        ## and a matrix nor with all "independent contributions" of the unique
        ## SKU-combinations. More details in our article.
        log_results == true ? print("\n  starting chi-square tests.") : nothing
        dep = HYOPTHESISTEST_SPARSE(Q, I, J, sig, ordered_skus, min_effect)
        ## Independence gate: compare the dependency premium carried by
        ## well-supported flagged pairs against the mean of the same
        ## statistic on ten margin-preserving randomizations of the
        ## transactions. Since each randomization preserves all order sizes
        ## and SKU frequencies, the total coappearance mass is identical and
        ## the masses compare directly; averaging stabilizes the chance
        ## level, whose single-draw variance is large when few pairs carry
        ## support. If the observed mass is within gate_ratio times its
        ## chance level, the flagged dependencies are statistically
        ## indistinguishable from noise and the allocation proceeds on
        ## expected coappearances alone.
        if gate_ratio > 0.0
            log_results == true ? print("\n  evaluating independence gate.") : nothing
            null_mass = 0.0
            for r in 1:10
                Q0 = COAPPEARENCE(CURVEBALL(trans; seed = 4242 + r), sku_weight)
                ordered_skus0 = [Q0[i, i] for i in axes(Q0, 1)]
                CLEANPRINCIPLE!(Q0)
                dep0 = HYOPTHESISTEST_SPARSE(Q0, I, J, sig, ordered_skus0, min_effect)
                null_mass += SUPPORTEDDEPMASS(dep0, Q0)
            end
            if SUPPORTEDDEPMASS(dep, Q) <= gate_ratio * null_mass / 10
                dep = sparse(Int32[], Int32[], Float32[], I, I)
                gate_fired = true
            end
        end
        ## Allocate every SKU once on capacities that hold the buffer free in every warehouse
        ## but the largest, then fill the capacity left. If the single copies cannot be
        ## packed, their target grows by the largest storage requirement.
        X = zeros(Bool, I, size(capacity, 1))
        attempt = 0
        packing_error = ErrorException("CHI: Cannot allocate every SKU once.")
        while true
            target = ceil(Int64, sum(sku_weight) + attempt * maximum(sku_weight))
            if attempt > 0 && target >= sum(capacity)
                throw(packing_error)
            end
            alloc_capacity = RESERVE(capacity, target)
            if attempt > 0
                print("\n  CHI: single copies repeated with target ", target, " (attempt ", attempt, ").")
            end
            p = sortperm(alloc_capacity; rev = true)
            packed = false
            try
                X, cap_left, nor_order = SINGLECOPIES(
                    dep, Q, alloc_capacity[p], sku_weight, ordered_skus, log_results
                )
                log_results == true ? print("\n  starting allocation to leftover space.") : nothing
                cap_left .+= (capacity .- alloc_capacity)[p]
                FILLRESERVE!(X, Q, cap_left, sku_weight, ordered_skus, nor_order, !gate_fired)
                X = X[:, invperm(p)]
                ## the tolerance absorbs rounding in the summed storage requirements
                packed = all(any(X; dims = 2)) && all(vec(X' * sku_weight) .<= capacity .* (1 + 1e-9))
            catch e
                if !(e isa ErrorException && occursin(r"^(SWAPREPAIR!|Cannot allocate SKU)", e.msg))
                    rethrow()
                end
                packing_error = e
            end
            packed && break
            attempt += 1
        end
        ## Free up RAM for last stage
        dep = nothing
        GC.gc()
        ls = 0
        if max_ls > 0 && !gate_fired
            if all(y->y==sku_weight[1], sku_weight)
                log_results == true ? print("\n  starting local search.") : nothing
                ls = LOCALSEARCHCHI!(
                    trans, X, Q, capacity, log_results, ls, max_ls, ls_neighborhood
                )
            else
                @warn "Local search requires uniform SKU weights — skipping. Results reflect CHI only."
            end
        end
        # Check whether all SKUs are allocated
        if any(y->y < 1, sum(X; dims = 2))
            error("CHI: Not all SKUs are allocated!")
        end
    end
    ## return the resulting allocation matrix
    return X, ls, gate_fired
end

# Allocation of every SKU to exactly one warehouse; capacity sorted decreasingly. Returns the
# allocation, the capacity left per warehouse and the order of the SKUs by their independent
# coappearances per storage unit (ties by order frequency).
function SINGLECOPIES(
    dep::SparseMatrixCSC{<:Real},
    Q::SparseMatrixCSC{<:Real},
    capacity::Vector{Int64},
    sku_weight::Vector{<:Real},
    ordered_skus::Vector{<:Real},
    log_results::Bool,
)
    I = size(Q, 1)
    ## Create the sku-warehouse allocation matrix
    X = zeros(Bool, I, size(capacity, 1))
    ## Determine the sum of the coappearances for each SKU on the base of
    ## the matrices dep and nor
    sum_dep = vec(sum(dep; dims = 2)) ./ sku_weight
    sum_nor = (vec(sum(Q; dims = 2)) .- vec(sum(dep; dims = 2))) ./ sku_weight
    ## Copy the capacity to have an array that keeps the left-over capacity
    cap_left::Vector{Float64} = copy(capacity)
    ## Create arrays that save the current dependencies for all unallocated
    ## SKUs to all warehouses and a vector that holds all unallocated skus
    state_dep = zeros(Float64, size(X, 1), length(cap_left))
    state_nor = zeros(Float64, size(X, 1), length(cap_left))
    allocated = zeros(Bool, size(X, 1))
    ## Calculate the weight of each warehouse. It shows us the density of
    ## the independent coappearances in each warehouse if we were to allocate
    ## all SKUs simply according to the highest independent coappearances.
    weight = WHWEIGHT(capacity, sum_nor, sku_weight)
    ## In the heuristic we compare the split-delivery minimising potential
    ## of SKU allocations to maximise dependent coappearances (SKUs with
    ## positive dependencies are stored in the same warehouse) with the
    ## maximisation of independent coappearances (SKUs with high independent
    ## coappearances are stored in warehouses with high capacities).
    ## Start the allocation:
    ## If all SKUs are allocated once, continue after the while loop, else
    log_results == true ? print("\n  starting unique allocation of each sku.") : nothing
    ## Allocate the SKUs with coapperances
    n_allocated = Ref(0)
    nor_order = sortperm(
        1:I; by = i -> (sum_nor[i], Float64(ordered_skus[i])), rev = true
    )
    nor_pos = Ref(1)
    while n_allocated[] < I
        ## select the SKU i with the highest coappearance not being allocated
        ## to the warehouses yet. In addition, select the first sorted
        ## warehouse k which has unused capacity left.
        ## Then, check for each warehouse with free space whether SKU i has
        ## significant dependencies to other already allocated SKUs.
        ## If dependencies exist check whether they are expected to dominate.
        ## If so, allocate SKU i to to the warehouse with the dependencies
        ## to complement the positive dependent SKUs in the warehouse.
        ## Otherwise allocate SKU i to warehouse k as the independent
        ## coappearances are expected to dominate.
        ## If no dependencies to yet allocated SKUs exist, select the warehouse
        ## with the highest unused capacity. If the dependet coappearances are
        ## expected to to dominate later, we allocate SKU i to the warehouse
        ## with the highest unused capacity. This maximises the likelihood to
        ## allocate all significant SKUs from a SKU-cluster into one warehouse.
        ## Otherwise we allocate it to warehouse k to maximise the independent
        ## coappearances in the allocation.
        i, k = SELECTIK(
            sum_dep,
            sum_nor,
            weight,
            cap_left,
            state_dep,
            allocated,
            sku_weight,
            nor_order,
            nor_pos,
        )
        ALLOCATEONE!(
            X,
            dep,
            Q,
            sum_dep,
            sum_nor,
            state_dep,
            state_nor,
            cap_left,
            allocated,
            sku_weight,
            i,
            k,
            n_allocated,
        )
        ## Check for all unallocated SKUs whether they have positive
        ## dependencies to the SKUs in the warehouse k the last SKU was allocated
        ## to. If so, check whether the dependencies are expected to dominate the
        ## independent coapperances. If yes, allocate the corresponding SKUs to
        ## the warehouse k.
        ADDDEPENDENT!(
            X,
            Q,
            cap_left,
            k,
            dep,
            sum_dep,
            sum_nor,
            state_dep,
            state_nor,
            allocated,
            sku_weight,
            n_allocated,
        )
        FILLLAST!(X, cap_left, allocated, sku_weight, n_allocated)
    end
    return X, cap_left, nor_order
end
