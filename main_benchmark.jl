function BENCHMARK(
    capacity_benchmark::Array{Int64,2},
    skus_benchmark::Vector{Int64},
    diff_benchmark::Vector{Float64},
    buff_benchmark::Vector{Float64},
    start::DataFrame,
    order_sets::Vector{Int64},
    trials::Int64,
    stagnant::Int64,
    strategy::Int64,
    klinkstatus::Bool,
    abort::Int64,
    iterations::Int64,
    show_opt::Bool,
    cpu_cores::Int64,
    allowed_gap::Float64,
    max_nodes::Int64,
    sig_levels::Vector{Float64},
    max_ls::Int64,
    chistatus::Bool,
    min_effect::Float64,
    ls_neighborhood::Float64,
    apply_ls::Bool,
    weight_mode::Symbol,
    weight_range::UnitRange{Int64},
    max_iih_iterations::Int64,
    epsilon_iih::Float64,
    benchiterations::Int64,
    train_test::Float64,
    dependency::String;
    gate_ratio::Float64 = 2.0,
    only_pairs::Union{Nothing,Set{Tuple{Int,Int}}} = nothing,
)

    # Open log file for this benchmark run
    logfile = "results/benchmark_$(experiment)_$(dependency).log"
    logio = open(logfile, "a")
    println(logio, "\n", "="^60)
    println(logio, "Benchmark started at ", Dates.now())
    println(logio, "Experiment: ", experiment, " | Dependency: ", dependency)
    println(logio, "="^60)
    flush(logio)

    function log_failure(
        io::IO, mode::String, e, bt, benchnr, skus, capacity, eff_capacity, orders, wmode
    )
        println(io, "\n[", Dates.now(), "] FAILURE: ", mode)
        println(
            io,
            "  benchiter=",
            benchnr,
            " skus=",
            skus,
            " wareh=",
            length(capacity),
            " capacity=",
            capacity,
            " effective_capacity=",
            eff_capacity,
            " orders=",
            orders,
            " weight_mode=",
            wmode,
        )
        println(io, "  Error: ", e)
        for line in split(sprint(showerror, e, bt), '\n')[1:min(end, 20)]
            println(io, "  ", line)
        end
        flush(io)
    end

    # Start the benchmark of the data set
    ## Create dataframes for the export of the results
    benchmark = DataFrame(;
        dependency = String[],
        skus = Int64[],
        wareh = Int64[],
        diff = Float64[],
        buffer = Float64[],
        weight_mode = String[],
        mode = String[],
        benchiter = Int64[],
        orders = Int64[],
        train_test = Float64[],
        parcel_train = Int64[],
        parcel_test = Int64[],
        flexibility = Float64[],
        duration = Float64[],
        cap_used = Int64[],
        local_search = Int64[],
        gap = Float64[],
    )

    ### Preallocate a transactional data set container
    trans = spzeros(Bool, 0, 0)

    # Repeat the benchmark benchiter times
    for benchnr in 1:benchiterations
        print("\nStarting benchmark iteration ", benchnr, " of ", benchiterations)

        # Iterate all capacity constellations
        for a in axes(capacity_benchmark, 1)
            # Restrict to selected (constellation, order_set) pairs; the
            # constellation index a stays the ORIGINAL row index so the
            # deterministic instance seeding matches the full runs
            if only_pairs !== nothing && !any(os -> (a, os) in only_pairs, order_sets)
                continue
            end
            ## Load the capacity of each individual run
            ### Note: It has to be sorted starting with the largest capacity
            capacity = Array{Int64,1}(undef, count(x -> x > 0, capacity_benchmark[a, :]))
            for b in axes(capacity, 1)
                capacity[b] = capacity_benchmark[a, b]
            end
            # Skip zero-buffer constellations for non-uniform weights
            if weight_mode != :uniform && buff_benchmark[a] <= 0.0
                @warn "Skipping constellation $a (zero buffer with non-uniform weights)"
                continue
            end

            print(
                "\n Benchmark ",
                benchnr,
                " of ",
                benchiterations,
                "/ Constellation: ",
                a,
                " of ",
                size(capacity_benchmark, 1),
                "\n Capacity: ",
                capacity,
            )

            ## Create all possible capacity combinations for the parcel Benchmark
            combination = COMBINEWAREHOUSES(capacity)

            ## Iterate over different number of orders
            for order_set in order_sets
                orders = order_set * skus_benchmark[a]
                if only_pairs !== nothing
                    (a, order_set) in only_pairs || continue
                    # force fresh generation: the reuse check below compares
                    # dimensions only, which is unsafe once pairs are skipped
                    trans = spzeros(Bool, 0, 0)
                end

                ## Generate artificial random transactions without dependencies if there is no transactional dataset
                if size(trans, 2) == skus_benchmark[a] && size(trans, 1) == orders
                    print("\n Reused transactions from previous run.")
                else
                    print("\n Starting generation of transactions.")
                    # Deterministic instance seeding: each scenario (scale,
                    # dependency, order density, constellation, repetition)
                    # regenerates the identical transactions on any machine,
                    # backing the reproducibility claim in the paper
                    Random.seed!(
                        skus_benchmark[a] + 31 * order_set + 1_009 * a +
                        65_537 * sum(Int, codeunits(dependency)) + 7 * benchnr,
                    )
                    max_gs = ceil(
                        Int64, max(skus_benchmark[a] / group_size_scaling, group_size_min)
                    )
                    mean_gs = ceil(Int64, max(max_gs / mean_group_divisor, mean_group_min))
                    time = @elapsed trans, C, group_sizes_gen = RANDOMTRANS(
                        skus_benchmark[a],
                        orders,
                        mean_order_size,
                        min_order_size,
                        nbd_dispersion,
                        sku_frequency_mode,
                        zipf_exponent,
                        max_gs,
                        mean_gs,
                        ratio_strong,
                        ratio_medium,
                        dep_strength_strong,
                        dep_strength_medium,
                        group_link,
                        one_direction,
                        multi_relatio,
                        dep_activation_prob,
                    )
                    print(
                        "\n Transactions generated after ",
                        round(time; digits = 3),
                        " seconds.",
                    )
                    if skus_benchmark[a] <= 10000
                        display(histogram(sum(trans; dims = 2)))
                        display(histogram(vec(sum(trans; dims = 1))))
                    end
                end

                #  Split the data into training and test data
                if train_test > 0.00
                    cut = round(Int64, size(trans, 1) * train_test)
                    trans_train = trans[1:cut, :]
                    trans_test = trans[(cut + 1):size(trans, 1), :]
                else
                    trans_train = trans_test = trans
                end
                print("\n Number of transactions for training ", size(trans_train, 1), ".")
                print("\n Number of transactions for validation ", size(trans_test, 1), ".")

                # SKU weight
                if weight_mode == :frequency
                    sku_weight = vec(Int64.(sum(trans_train; dims = 1)))
                    sku_weight = max.(sku_weight, 1)  # ensure minimum weight of 1
                elseif weight_mode == :random
                    rng_w = MersenneTwister(benchnr * 1000 + a * 100 + order_set)
                    sku_weight = rand(rng_w, weight_range, size(trans_train, 2))
                else
                    sku_weight = ones(Int64, size(trans_train, 2))
                end

                # Scale capacity for non-uniform weights
                if weight_mode != :uniform
                    avg_weight = sum(sku_weight) / length(sku_weight)
                    # Scale factor must ensure largest warehouse >= heaviest SKU
                    scale = max(avg_weight, maximum(sku_weight) / maximum(capacity))
                    effective_capacity = round.(Int64, capacity .* scale)
                    # Ensure total capacity >= total weight
                    while sum(effective_capacity) < sum(sku_weight)
                        effective_capacity[argmax(effective_capacity)] += 1
                    end
                    effective_combination = COMBINEWAREHOUSES(effective_capacity)
                    print(
                        "\n Weight mode: ",
                        weight_mode,
                        " / avg weight: ",
                        round(avg_weight; digits = 2),
                        " / effective capacity: ",
                        effective_capacity,
                    )
                else
                    effective_capacity = capacity
                    effective_combination = combination
                end

                # Compute coappearance matrix for general local search
                can_ls = apply_ls && all(y->y==sku_weight[1], sku_weight)
                time_Q_ls = 0.0
                if can_ls
                    time_Q_ls = @elapsed begin
                        Q_ls = COAPPEARENCE(trans_train, sku_weight)
                        CLEANPRINCIPLE!(Q_ls)
                    end
                end

                #  Start the heuristics and optmisations

                ## Start QMK heuristic with Gurobi as solver
                if start[1, :QMK] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, gap_optimisation = MQKP(
                            trans_train,
                            effective_capacity,
                            sku_weight,
                            abort,
                            "Gurobi",
                            show_opt,
                            cpu_cores,
                            allowed_gap,
                            max_nodes,
                            "QMK",
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n    QMK: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / gap: ",
                            round(gap_optimisation; digits = 6),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "QMK",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "QMK failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "QMK",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start QMK heuristic with Juniper as solver
                if start[1, :QMKJ] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, gap_optimisation = MQKP(
                            trans_train,
                            effective_capacity,
                            sku_weight,
                            abort,
                            "Juniper",
                            show_opt,
                            cpu_cores,
                            allowed_gap,
                            max_nodes,
                            "QMK",
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n   QMKJ: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / gap: ",
                            round(gap_optimisation; digits = 6),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "QMKJ",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "QMKJ failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "QMKJ",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start QMK heuristic with SCIP as solver
                if start[1, :QMKS] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, gap_optimisation = MQKP(
                            trans_train,
                            effective_capacity,
                            sku_weight,
                            abort,
                            "scip",
                            show_opt,
                            cpu_cores,
                            allowed_gap,
                            max_nodes,
                            "QMK",
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n   QMKS: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / gap: ",
                            round(gap_optimisation; digits = 6),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "QMKS",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "QMKS failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "QMKS",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start chi square heuristic without local search
                if start[1, :CHI] == 1
                    try
                        for sig in sig_levels
                            sleep(0.01)
                            GC.gc()
                            time_benchmark = @elapsed W, ls = CHISQUAREHEUR(
                                trans_train,
                                effective_capacity,
                                sig,
                                0,
                                sku_weight,
                                chistatus,
                                min_effect,
                                ls_neighborhood,
                                gate_ratio,
                            )
                            parcels_benchmark = PARCELSSEND(
                                trans_test, W, effective_capacity, effective_combination
                            )
                            parcels_train = PARCELSSEND(
                                trans_train, W, effective_capacity, effective_combination
                            )
                            flex_benchmark = FLEXIBILITY(trans_test, W)
                            print(
                                "\n   CHI: parcels test data: ",
                                parcels_benchmark,
                                " / parcels training data: ",
                                parcels_train,
                                " / flex: ",
                                flex_benchmark,
                                " / time: ",
                                round(time_benchmark; digits = 3),
                                " / local search: ",
                                ls,
                                " / sig: ",
                                sig,
                                " / warehouse: ",
                                sum(W; dims = 1),
                            )

                            push!(
                                benchmark,
                                (
                                    dependency = dependency,
                                    skus = skus_benchmark[a],
                                    wareh = length(capacity),
                                    diff = diff_benchmark[a],
                                    buffer = buff_benchmark[a],
                                    weight_mode = string(weight_mode),
                                    mode = "CHI_$sig",
                                    benchiter = benchnr,
                                    orders = size(trans, 1),
                                    train_test = train_test,
                                    parcel_train = parcels_train,
                                    parcel_test = parcels_benchmark,
                                    flexibility = flex_benchmark,
                                    duration = time_benchmark,
                                    cap_used = sum(W),
                                    local_search = ls,
                                    gap = 0,
                                ),
                            )
                            if can_ls
                                RUN_LS!(
                                    benchmark,
                                    W,
                                    trans_train,
                                    trans_test,
                                    Q_ls,
                                    effective_capacity,
                                    effective_combination,
                                    max_ls,
                                    "CHI_$sig",
                                    dependency,
                                    skus_benchmark[a],
                                    length(capacity),
                                    diff_benchmark[a],
                                    buff_benchmark[a],
                                    string(weight_mode),
                                    benchnr,
                                    size(trans, 1),
                                    train_test,
                                    time_Q_ls,
                                    time_benchmark,
                                )
                            end
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "CHI failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "CHI",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start our reproduction of the  K-LINKS heuristic by
                ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
                if start[1, :KL] == 1 &&
                    all(w->w==sku_weight[1], sku_weight) &&
                    sum(effective_capacity) == sum(sku_weight)
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, ls = KLINKS(
                            trans_train,
                            effective_capacity,
                            trials,
                            stagnant,
                            strategy,
                            abort,
                            klinkstatus,
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n     KL: parcels train data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / local search: ",
                            ls,
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "KL",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = ls,
                                gap = 0,
                            ),
                        )
                        if can_ls
                            RUN_LS!(
                                benchmark,
                                W,
                                trans_train,
                                trans_test,
                                Q_ls,
                                effective_capacity,
                                effective_combination,
                                max_ls,
                                "KL",
                                dependency,
                                skus_benchmark[a],
                                length(capacity),
                                diff_benchmark[a],
                                buff_benchmark[a],
                                string(weight_mode),
                                benchnr,
                                size(trans, 1),
                                train_test,
                                time_Q_ls,
                                time_benchmark,
                            )
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "KL failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "KL",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start our reproduction of the  K-LINKS optimization with SBB by
                ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
                if start[1, :KLQ] == 1 &&
                    all(w->w==sku_weight[1], sku_weight) &&
                    sum(effective_capacity) == sum(sku_weight)
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, gap_optimisation = MQKP(
                            trans_train,
                            effective_capacity,
                            sku_weight,
                            abort,
                            "Gurobi",
                            show_opt,
                            cpu_cores,
                            allowed_gap,
                            max_nodes,
                            "QMK",
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n    KLQ: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / gap: ",
                            round(gap_optimisation; digits = 6),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "KLQ",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "KLQ failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "KLQ",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start our reproduction of the greedy orders heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :GO] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W = GREEDYORDERS(
                            trans_train, effective_capacity, sku_weight
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n     GO: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "GO",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = 0,
                            ),
                        )
                        if can_ls
                            RUN_LS!(
                                benchmark,
                                W,
                                trans_train,
                                trans_test,
                                Q_ls,
                                effective_capacity,
                                effective_combination,
                                max_ls,
                                "GO",
                                dependency,
                                skus_benchmark[a],
                                length(capacity),
                                diff_benchmark[a],
                                buff_benchmark[a],
                                string(weight_mode),
                                benchnr,
                                size(trans, 1),
                                train_test,
                                time_Q_ls,
                                time_benchmark,
                            )
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "GO failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "GO",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start our reproduction of the greedy pairs heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :GP] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W = GREEDYPAIRS(
                            trans_train, effective_capacity, sku_weight
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n     GP: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "GP",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = 0,
                            ),
                        )
                        if can_ls
                            RUN_LS!(
                                benchmark,
                                W,
                                trans_train,
                                trans_test,
                                Q_ls,
                                effective_capacity,
                                effective_combination,
                                max_ls,
                                "GP",
                                dependency,
                                skus_benchmark[a],
                                length(capacity),
                                diff_benchmark[a],
                                buff_benchmark[a],
                                string(weight_mode),
                                benchnr,
                                size(trans, 1),
                                train_test,
                                time_Q_ls,
                                time_benchmark,
                            )
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "GP failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "GP",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start our reproduction of the greedy seeds heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :GS] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W = GREEDYSEEDS(
                            trans_train, effective_capacity, sku_weight
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n     GS: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "GS",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = 0,
                            ),
                        )
                        if can_ls
                            RUN_LS!(
                                benchmark,
                                W,
                                trans_train,
                                trans_test,
                                Q_ls,
                                effective_capacity,
                                effective_combination,
                                max_ls,
                                "GS",
                                dependency,
                                skus_benchmark[a],
                                length(capacity),
                                diff_benchmark[a],
                                buff_benchmark[a],
                                string(weight_mode),
                                benchnr,
                                size(trans, 1),
                                train_test,
                                time_Q_ls,
                                time_benchmark,
                            )
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "GS failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "GS",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start our reproduction of the  bestselling heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :BS] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W = BESTSELLING(
                            trans_train, effective_capacity, sku_weight
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n     BS: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "BS",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = 0,
                            ),
                        )
                        if can_ls
                            RUN_LS!(
                                benchmark,
                                W,
                                trans_train,
                                trans_test,
                                Q_ls,
                                effective_capacity,
                                effective_combination,
                                max_ls,
                                "BS",
                                dependency,
                                skus_benchmark[a],
                                length(capacity),
                                diff_benchmark[a],
                                buff_benchmark[a],
                                string(weight_mode),
                                benchnr,
                                size(trans, 1),
                                train_test,
                                time_Q_ls,
                                time_benchmark,
                            )
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "BS failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "BS",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start the extended MCI heuristic for D warehouses by
                ## Lin et al. (2025) https://doi.org/10.1111/poms.14114
                if start[1, :EMCI] == 1
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W = EMCIALLOC(
                            trans_train, effective_capacity, sku_weight
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n   EMCI: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "EMCI",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = 0,
                            ),
                        )
                        if can_ls
                            RUN_LS!(
                                benchmark,
                                W,
                                trans_train,
                                trans_test,
                                Q_ls,
                                effective_capacity,
                                effective_combination,
                                max_ls,
                                "EMCI",
                                dependency,
                                skus_benchmark[a],
                                length(capacity),
                                diff_benchmark[a],
                                buff_benchmark[a],
                                string(weight_mode),
                                benchnr,
                                size(trans, 1),
                                train_test,
                                time_Q_ls,
                                time_benchmark,
                            )
                        end
                    catch e
                        bt = catch_backtrace()
                        @warn "EMCI failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "EMCI",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start the iterative improvement heuristic with Gurobi by
                ## Lin et al. (2025) https://doi.org/10.1111/poms.14114
                if start[1, :IIH] == 1 &&
                    all(w->w==sku_weight[1], sku_weight) &&
                    length(effective_capacity) == 2 &&
                    sum(effective_capacity) > sum(sku_weight)
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, gap_optimisation, ls = IIH(
                            trans_train,
                            effective_capacity,
                            sku_weight,
                            abort,
                            "Gurobi",
                            show_opt,
                            cpu_cores,
                            allowed_gap,
                            max_nodes,
                            max_iih_iterations,
                            epsilon_iih,
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n    IIH: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / gap: ",
                            round(gap_optimisation; digits = 6),
                            " / iterations: ",
                            ls,
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "IIH",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = ls,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "IIH failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "IIH",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start the iterative improvement heuristic with SCIP by
                ## Lin et al. (2025) https://doi.org/10.1111/poms.14114
                if start[1, :IIHS] == 1 &&
                    all(w->w==sku_weight[1], sku_weight) &&
                    length(effective_capacity) == 2 &&
                    sum(effective_capacity) > sum(sku_weight)
                    try
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, gap_optimisation, ls = IIH(
                            trans_train,
                            effective_capacity,
                            sku_weight,
                            abort,
                            "scip",
                            show_opt,
                            cpu_cores,
                            allowed_gap,
                            max_nodes,
                            max_iih_iterations,
                            epsilon_iih,
                        )
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n   IIHS: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / gap: ",
                            round(gap_optimisation; digits = 6),
                            " / iterations: ",
                            ls,
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "IIHS",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = ls,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "IIHS failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "IIHS",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Start the search for optimal solution with the solver CPLEX
                ## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
                ## FULLOPTUEQ if SKUs can be allocated multiple times
                if start[1, :OPT] == 1 && all(w->w==sku_weight[1], sku_weight)
                    try
                        sleep(0.01)
                        GC.gc()
                        if sum(effective_capacity) == size(trans, 2)
                            time_benchmark = @elapsed W, gap_optimisation, popt = FULLOPTEQ(
                                trans_train,
                                effective_capacity,
                                abort,
                                show_opt,
                                cpu_cores,
                                allowed_gap,
                                max_nodes,
                            )
                        else
                            time_benchmark = @elapsed W, gap_optimisation, popt = FULLOPTUEQ(
                                trans_train,
                                effective_capacity,
                                abort,
                                show_opt,
                                cpu_cores,
                                allowed_gap,
                                max_nodes,
                            )
                        end
                        parcels_benchmark = PARCELSSEND(
                            trans_test, W, effective_capacity, effective_combination
                        )
                        parcels_train = PARCELSSEND(
                            trans_train, W, effective_capacity, effective_combination
                        )
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print(
                            "\n    OPT: parcels test data: ",
                            parcels_benchmark,
                            " / parcels training data: ",
                            parcels_train,
                            " / flex: ",
                            flex_benchmark,
                            " / time: ",
                            round(time_benchmark; digits = 3),
                            " / warehouse: ",
                            sum(W; dims = 1),
                        )

                        push!(
                            benchmark,
                            (
                                dependency = dependency,
                                skus = skus_benchmark[a],
                                wareh = length(capacity),
                                diff = diff_benchmark[a],
                                buffer = buff_benchmark[a],
                                weight_mode = string(weight_mode),
                                mode = "OPT",
                                benchiter = benchnr,
                                orders = size(trans, 1),
                                train_test = train_test,
                                parcel_train = parcels_train,
                                parcel_test = parcels_benchmark,
                                flexibility = flex_benchmark,
                                duration = time_benchmark,
                                cap_used = sum(W),
                                local_search = 0,
                                gap = gap_optimisation,
                            ),
                        )
                    catch e
                        bt = catch_backtrace()
                        @warn "OPT failed" exception=(e, bt)
                        log_failure(
                            logio,
                            "OPT",
                            e,
                            bt,
                            benchnr,
                            skus_benchmark[a],
                            capacity,
                            effective_capacity,
                            orders,
                            weight_mode,
                        )
                    end
                end

                ## Benchmark the random allocation of SKUs
                try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed parcels_benchmark, flex_benchmark = RANDOMBENCH(
                        trans_test,
                        effective_capacity,
                        iterations,
                        sku_weight,
                        effective_combination,
                    )
                    parcels_train, _ = RANDOMBENCH(
                        trans_train,
                        effective_capacity,
                        iterations,
                        sku_weight,
                        effective_combination,
                    )
                    print(
                        "\n    RND: parcels test data: ",
                        parcels_benchmark,
                        " / parcels training data: ",
                        parcels_train,
                        " / flex: ",
                        flex_benchmark,
                        " / time: ",
                        round(time_benchmark; digits = 3),
                        " / benchmarks: ",
                        iterations,
                        " / warehouse: ",
                        sum(W; dims = 1),
                        "\n",
                    )

                    push!(
                        benchmark,
                        (
                            dependency = dependency,
                            skus = skus_benchmark[a],
                            wareh = length(capacity),
                            diff = diff_benchmark[a],
                            buffer = buff_benchmark[a],
                            weight_mode = string(weight_mode),
                            mode = "RND",
                            benchiter = benchnr,
                            orders = size(trans, 1),
                            train_test = train_test,
                            parcel_train = parcels_train,
                            parcel_test = parcels_benchmark,
                            flexibility = flex_benchmark,
                            duration = time_benchmark,
                            cap_used = sum(W),
                            local_search = 0,
                            gap = 0,
                        ),
                    )
                catch e
                    bt = catch_backtrace()
                    @warn "RND failed" exception=(e, bt)
                    log_failure(
                        logio,
                        "RND",
                        e,
                        bt,
                        benchnr,
                        skus_benchmark[a],
                        capacity,
                        effective_capacity,
                        orders,
                        weight_mode,
                    )
                end

                # Free training/test data before next iteration
                trans_train = nothing
                trans_test = nothing
                sku_weight = nothing
                W = nothing
                can_ls && (Q_ls = nothing)
                GC.gc()

                # Export the results after each stage
                wm_suffix = weight_mode == :uniform ? "" : "_$(weight_mode)"
                CSV.write(
                    "results/$(experiment)_benchmark_$(dependency)$(wm_suffix).csv",
                    benchmark,
                )
            end
        end
    end

    print("\n### Finished Benchmark ###")

    println(logio, "\n", "="^60)
    println(logio, "Benchmark finished at ", Dates.now())
    println(logio, "="^60)
    close(logio)

    return benchmark::DataFrame
end
