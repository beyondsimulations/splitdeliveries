# Case-study benchmark functions used by apply_cbench.jl.

## PARCELSSEND_WEIGHT: split evaluation that additionally reports how many
## orders are dispatched (fully or partially) from each warehouse.
## Each order is served by the smallest warehouse combination that covers it;
## prefer_last flips the tie-break among equally small combinations, so calling
## the function with true and false yields the minimum and maximum number of
## dispatches per warehouse (the split count itself is identical).
## Unlike PARCELSSEND, SKUs without any allocation are tolerated: an order
## containing a SKU that is stocked nowhere (the warehouse simulation
## evaluates daily snapshots in which seasonal SKUs are missing) is dispatched
## from the full warehouse set, preserving the identity
## sum(dispatches) = orders + splits throughout.
## Returns (splits, dispatches): splits = parcels beyond one per order;
## dispatches[k] = orders served (at least partly) from warehouse k.
function PARCELSSEND_WEIGHT(
    trans::SparseMatrixCSC{Bool,Int64},
    X::Array{Bool,2},
    capacity::Vector{<:Real},
    combination::Array{Array{Array{Int64,1},1},1},
    prefer_last::Bool,
)

    # SKU coverage of each warehouse combination
    ncomb = size(combination, 1)
    warehouses = [Int64[] for _ in 1:ncomb]
    warehouse_combination = zeros(Int64, size(X, 1), ncomb)
    for i in 1:ncomb
        for j in axes(combination[i], 1)
            for k in combination[i][j]
                push!(warehouses[i], k)
                for s in axes(X, 1)
                    if X[s, k]
                        warehouse_combination[s, i] = 1
                    end
                end
            end
        end
    end
    sizes = [length(warehouses[i]) for i in 1:ncomb]

    # Scan order: smallest combinations first, tie-break flipped by prefer_last.
    # combinations() may emit the empty set; exclude it, otherwise orders
    # without any covered SKU would be served by no warehouse at all.
    scan = sort([i for i in 1:ncomb if sizes[i] > 0]; by = i -> (sizes[i], prefer_last ? -i : i))

    fulfilled = Matrix(trans * dropzeros(sparse(warehouse_combination)))
    need = vec(sum(trans; dims = 2))
    dispatches = zeros(Int64, size(X, 2))
    parcel = 0
    for o in axes(fulfilled, 1)
        served = false
        for i in scan
            if fulfilled[o, i] == need[o]
                parcel += sizes[i]
                for k in warehouses[i]
                    dispatches[k] += 1
                end
                served = true
                break
            end
        end
        if !served
            # the order contains SKUs stocked nowhere in this allocation:
            # dispatch from the full warehouse set
            full = scan[end]
            parcel += sizes[full]
            for k in warehouses[full]
                dispatches[k] += 1
            end
        end
    end
    parcel -= size(trans, 1)
    return parcel, dispatches
end

## benchmark!: run all algorithms enabled in `start` on one test day of the
## case study and append one result row per algorithm (plus RND) to `benchmark`.
## Allocations persist in X (categorybrands x warehouses x algorithm) between
## days; the heuristics are re-run only when optimization_cycle is true (or on
## the first day), mirroring the weekly re-optimisation of the case study.
## `threshold` is kept for signature compatibility with an earlier
## allocation-smoothing experiment and is not used.
function benchmark!(
    start::DataFrame,
    X::Array{Bool,3},
    benchmark::DataFrame,
    categorybrands::DataFrame,
    datasource::String,
    training_days::Int64,
    date,
    train_orders::SparseMatrixCSC{Bool,Int64},
    test_orders::SparseMatrixCSC{Bool,Int64},
    capacity::Vector{Int64},
    abort::Int64,
    show_opt::Bool,
    cpu_cores::Int64,
    allowed_gap::Float64,
    max_nodes::Int64,
    sig::Float64,
    chistatus::Bool,
    max_ls::Int64,
    trials::Int64,
    stagnant::Int64,
    strategy::Int64,
    klinkstatus::Bool,
    sku_weight::Vector{<:Real},
    optimization_cycle::Bool,
    training_intervall::Int64,
    threshold::Int64;
    min_effect::Float64 = 0.0,
    ls_neighborhood::Float64 = 1.0,
)
    if length(capacity) != 2
        error("benchmark!: the case study expects exactly two warehouses.")
    end
    dict_algorithm = Dict(String(name) => i for (i, name) in enumerate(names(start)))
    combination = COMBINEWAREHOUSES(capacity)

    # The heuristics sort capacities in decreasing order internally; pass them
    # sorted and map the allocation columns back to the input warehouse order
    perm = sortperm(capacity; rev = true)
    cap_sorted = capacity[perm]
    unsort = invperm(perm)

    diff = round((maximum(capacity) - minimum(capacity)) / sum(capacity); digits = 2)
    buffer = round(100 * (sum(capacity) - sum(sku_weight)) / sum(sku_weight); digits = 2)
    uniform = all(w -> w == sku_weight[1], sku_weight)

    enabled(alg) = alg in names(start) && start[1, alg] == 1

    ## Evaluate one allocation on the current train/test windows, append a row
    function record!(mode::String, W::Matrix{Bool}, duration, reoCB, reoSKU, ls, gap)
        parcels_train = PARCELSSEND(train_orders, W, capacity, combination)
        parcels_test, disp_high = PARCELSSEND_WEIGHT(
            test_orders, W, capacity, combination, true
        )
        parcels_test, disp_low = PARCELSSEND_WEIGHT(
            test_orders, W, capacity, combination, false
        )
        weight_app = floor.(Int64, vec(sum(W .* sku_weight; dims = 1)))
        print(
            "\n ", rpad(mode, 10),
            ": splits test: ", parcels_test,
            " / splits train: ", parcels_train,
            " / time: ", round(duration; digits = 3),
            " / reorders: ", reoSKU,
            " / warehouse: ", sum(W; dims = 1),
        )
        push!(
            benchmark,
            (
                data = datasource,
                skus = nrow(categorybrands),
                orders_train = size(train_orders, 1),
                orders_test = size(test_orders, 1),
                date = date,
                traindays = training_days,
                trainingintervall = training_intervall,
                capacity_47 = capacity[1],
                capacity_50 = capacity[2],
                min_dispatch_47 = disp_high[1],
                min_dispatch_50 = disp_low[2],
                max_dispatch_47 = disp_low[1],
                max_dispatch_50 = disp_high[2],
                weight_app_47 = weight_app[1],
                weight_app_50 = weight_app[2],
                diff = diff,
                buffer = buffer,
                mode = mode,
                parcel_train = parcels_train,
                parcel_test = parcels_test,
                reordercategorybrands = reoCB,
                reorderskus = float(reoSKU),
                duration = float(duration),
                cap_used = sum(W),
                local_search = ls,
                gap = float(gap),
            ),
        )
    end

    ## Run one algorithm respecting the re-optimisation cycle: `solve` returns
    ## (W, ls, gap) with W in the input warehouse order
    function cycle!(alg::String, mode::String, solve::Function)
        idx = dict_algorithm[alg]
        if optimization_cycle || sum(@view X[:, :, idx]) == 0
            sleep(0.01)
            GC.gc()
            local W, ls, gap
            duration = @elapsed W, ls, gap = solve()
            reoCB, reoSKU, _ = reorders(X, W, dict_algorithm, categorybrands, alg)
            record!(mode, W, duration, reoCB, reoSKU, ls, gap)
        else
            W = Matrix{Bool}(X[:, :, idx])
            record!(mode, W, 0.0, 0, 0.0, 0, 0.0)
        end
    end

    ## Quadratic multiple knapsack heuristic (Gurobi)
    if enabled("QMKO")
        try
            cycle!(
                "QMKO",
                "QMKO",
                () -> begin
                    W, gap = MQKP(
                        train_orders, cap_sorted, sku_weight, abort, "Gurobi",
                        show_opt, cpu_cores, allowed_gap, max_nodes, "QMK",
                    )
                    (W[:, unsort], 0, gap)
                end,
            )
        catch e
            @warn "QMKO failed" exception = (e, catch_backtrace())
        end
    end

    ## Quadratic multiple knapsack heuristic (Juniper)
    if enabled("QMK")
        try
            cycle!(
                "QMK",
                "QMK",
                () -> begin
                    W, gap = MQKP(
                        train_orders, cap_sorted, sku_weight, abort, "Juniper",
                        show_opt, cpu_cores, allowed_gap, max_nodes, "QMK",
                    )
                    (W[:, unsort], 0, gap)
                end,
            )
        catch e
            @warn "QMK failed" exception = (e, catch_backtrace())
        end
    end

    ## Chi-square heuristic without local search
    if enabled("CHIM")
        try
            cycle!(
                "CHIM",
                "CHIM_$sig",
                () -> begin
                    W, ls = CHISQUAREHEUR(
                        train_orders, cap_sorted, sig, 0, sku_weight,
                        chistatus, min_effect, ls_neighborhood,
                    )
                    (W[:, unsort], ls, 0.0)
                end,
            )
        catch e
            @warn "CHIM failed" exception = (e, catch_backtrace())
        end
    end

    ## Chi-square heuristic with local search (uniform SKU weights only)
    if enabled("CHI")
        if uniform
            try
                cycle!(
                    "CHI",
                    "CHI_$sig",
                    () -> begin
                        W, ls = CHISQUAREHEUR(
                            train_orders, cap_sorted, sig, max_ls, sku_weight,
                            chistatus, min_effect, ls_neighborhood,
                        )
                        (W[:, unsort], ls, 0.0)
                    end,
                )
            catch e
                @warn "CHI failed" exception = (e, catch_backtrace())
            end
        else
            @warn "CHI skipped: local search requires uniform SKU weights."
        end
    end

    ## K-LINKS heuristic (uniform weights and exactly matching capacity only)
    if enabled("KL")
        if uniform && sum(capacity) == sum(sku_weight)
            try
                cycle!(
                    "KL",
                    "KL",
                    () -> begin
                        W, ls = KLINKS(
                            train_orders, cap_sorted, trials, stagnant,
                            strategy, abort, klinkstatus,
                        )
                        (W[:, unsort], ls, 0.0)
                    end,
                )
            catch e
                @warn "KL failed" exception = (e, catch_backtrace())
            end
        else
            @warn "KL skipped: requires uniform SKU weights and capacity equal to the number of allocation objects."
        end
    end

    ## K-LINKS optimisation (uniform weights and exactly matching capacity only)
    if enabled("KLQ")
        if uniform && sum(capacity) == sum(sku_weight)
            try
                cycle!(
                    "KLQ",
                    "KLQ",
                    () -> begin
                        W, gap = MQKP(
                            train_orders, cap_sorted, sku_weight, abort, "Gurobi",
                            show_opt, cpu_cores, allowed_gap, max_nodes, "QMK",
                        )
                        (W[:, unsort], 0, gap)
                    end,
                )
            catch e
                @warn "KLQ failed" exception = (e, catch_backtrace())
            end
        else
            @warn "KLQ skipped: requires uniform SKU weights and capacity equal to the number of allocation objects."
        end
    end

    ## Greedy orders heuristic
    if enabled("GO")
        try
            cycle!(
                "GO",
                "GO",
                () -> (GREEDYORDERS(train_orders, cap_sorted, sku_weight)[:, unsort], 0, 0.0),
            )
        catch e
            @warn "GO failed" exception = (e, catch_backtrace())
        end
    end

    ## Greedy pairs heuristic
    if enabled("GP")
        try
            cycle!(
                "GP",
                "GP",
                () -> (GREEDYPAIRS(train_orders, cap_sorted, sku_weight)[:, unsort], 0, 0.0),
            )
        catch e
            @warn "GP failed" exception = (e, catch_backtrace())
        end
    end

    ## Greedy seeds heuristic
    if enabled("GS")
        try
            cycle!(
                "GS",
                "GS",
                () -> (GREEDYSEEDS(train_orders, cap_sorted, sku_weight)[:, unsort], 0, 0.0),
            )
        catch e
            @warn "GS failed" exception = (e, catch_backtrace())
        end
    end

    ## Bestsellers heuristic
    if enabled("BS")
        try
            cycle!(
                "BS",
                "BS",
                () -> (BESTSELLING(train_orders, cap_sorted, sku_weight)[:, unsort], 0, 0.0),
            )
        catch e
            @warn "BS failed" exception = (e, catch_backtrace())
        end
    end

    ## Extended MCI heuristic (optional column in `start`)
    if enabled("EMCI")
        try
            cycle!(
                "EMCI",
                "EMCI",
                () -> (EMCIALLOC(train_orders, cap_sorted, sku_weight)[:, unsort], 0, 0.0),
            )
        catch e
            @warn "EMCI failed" exception = (e, catch_backtrace())
        end
    end

    ## Iterative improvement heuristic (optional column; uniform weights, two
    ## warehouses, buffered capacity only)
    if enabled("IIH")
        if uniform && sum(capacity) > sum(sku_weight)
            try
                cycle!(
                    "IIH",
                    "IIH",
                    () -> begin
                        W, gap, ls = IIH(
                            train_orders, cap_sorted, sku_weight, abort, "Gurobi",
                            show_opt, cpu_cores, allowed_gap, max_nodes, 50, 1.0e-6,
                        )
                        (W[:, unsort], ls, gap)
                    end,
                )
            catch e
                @warn "IIH failed" exception = (e, catch_backtrace())
            end
        else
            @warn "IIH skipped: requires uniform SKU weights and buffered capacity."
        end
    end

    ## Exact optimisation (uniform SKU weights only)
    if enabled("OPT")
        if uniform
            try
                cycle!(
                    "OPT",
                    "OPT",
                    () -> begin
                        if sum(capacity) == size(train_orders, 2)
                            W, gap, _ = FULLOPTEQ(
                                train_orders, cap_sorted, abort, show_opt,
                                cpu_cores, allowed_gap, max_nodes,
                            )
                        else
                            W, gap, _ = FULLOPTUEQ(
                                train_orders, cap_sorted, abort, show_opt,
                                cpu_cores, allowed_gap, max_nodes,
                            )
                        end
                        (W[:, unsort], 0, gap)
                    end,
                )
            catch e
                @warn "OPT failed" exception = (e, catch_backtrace())
            end
        else
            @warn "OPT skipped: requires uniform SKU weights."
        end
    end

    ## Random allocation baseline (re-drawn every day)
    try
        sleep(0.01)
        GC.gc()
        local W_rnd
        # float capacities: the fractional case-study weights (mean SKUs per
        # categorybrand) would otherwise fail the integer capacity updates
        duration = @elapsed W_rnd = RANDOMALLOCMULTI(
            test_orders, float.(capacity), sku_weight
        )
        record!("RND", W_rnd, duration, 0, 0.0, 0, 0.0)
    catch e
        @warn "RND failed" exception = (e, catch_backtrace())
    end

    return benchmark
end
