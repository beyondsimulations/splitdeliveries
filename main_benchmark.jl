function BENCHMARK(capacity_benchmark::Array{Int64,2},
    skus_benchmark::Vector{Int64},
    diff_benchmark::Vector{Float64},
    buff_benchmark::Vector{Float64},
    start::DataFrame,
    order_sets::Vector{Int64},
    max_dependence::Float64,
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
    max_iih_iterations::Int64,
    epsilon_iih::Float64,
    benchiterations::Int64,
    train_test::Float64,
    dependency::String)

    # Open log file for this benchmark run
    logfile = "results/benchmark_$(dependency).log"
    logio = open(logfile, "a")
    println(logio, "\n", "="^60)
    println(logio, "Benchmark started at ", Dates.now())
    println(logio, "Dependency: ", dependency)
    println(logio, "="^60)
    flush(logio)

    function log_failure(io::IO, mode::String, e, bt, benchnr, skus, capacity, orders)
        println(io, "\n[", Dates.now(), "] FAILURE: ", mode)
        println(io, "  benchiter=", benchnr, " skus=", skus, " wareh=", length(capacity),
            " capacity=", capacity, " orders=", orders)
        println(io, "  Error: ", e)
        for line in split(sprint(showerror, e, bt), '\n')[1:min(end, 20)]
            println(io, "  ", line)
        end
        flush(io)
    end

    # Start the benchmark of the data set
    ## Create dataframes for the export of the results
    benchmark = DataFrame(dependency=String[],
        skus=Int64[],
        wareh=Int64[],
        diff=Float64[],
        buffer=Float64[],
        mode=String[],
        benchiter=Int64[],
        orders=Int64[],
        train_test=Float64[],
        parcel_train=Int64[],
        parcel_test=Int64[],
        flexibility=Float64[],
        duration=Float64[],
        cap_used=Int64[],
        local_search=Int64[],
        gap=Float64[])

    ### Preallocate a transactional data set container
    trans = spzeros(Bool, 0, 0)

    # Repeat the benchmark benchiter times
    for benchnr = 1:benchiterations
        print("\nStarting benchmark iteration ", benchnr, " of ", benchiterations)

        # Iterate all capacity constellations
        for a in axes(capacity_benchmark, 1)
            ## Load the capacity of each individual run
            ### Note: It has to be sorted starting with the largest capacity
            capacity = Array{Int64,1}(undef, count(x -> x > 0, capacity_benchmark[a, :]))
            for b in axes(capacity, 1)
                capacity[b] = capacity_benchmark[a, b]
            end
            print("\n Benchmark ", benchnr, " of ", benchiterations,
                "/ Constellation: ", a, " of ", size(capacity_benchmark, 1),
                "\n Capacity: ", capacity)

            ## Create all possible capacity combinations for the parcel Benchmark
            combination = COMBINEWAREHOUSES(capacity)

            ## Iterate over different number of orders
            for order_set in order_sets

                orders = order_set * skus_benchmark[a]

                ## Generate artificial random transactions without dependencies if there is no transactional dataset
                if size(trans, 2) == skus_benchmark[a] && size(trans, 1) == orders
                    print("\n Reused transactions from previous run.")
                else
                    print("\n Starting generation of transactions.")
                    time = @elapsed trans = RANDOMTRANS(skus_benchmark[a], orders, skus_in_order, sku_frequency,
                        ceil(Int64, max(skus_benchmark[a] / 20, 10)),
                        min_dependence, max_dependence,
                        group_link, ind_chance, one_direction,
                        multi_relatio)
                    print("\n Transactions generated after ", round(time, digits=3), " seconds.")
                    if skus_benchmark[a] <= 10000
                        display(histogram(sum(trans, dims=2)))
                        display(histogram(vec(sum(trans, dims=1))))
                    end
                end

                #  Split the data into training and test data
                if train_test > 0.00
                    cut = round(Int64, size(trans, 1) * train_test)
                    trans_train = trans[1:cut, :]
                    trans_test = trans[(cut+1):size(trans, 1), :]
                else
                    trans_train = trans_test = trans
                end
                print("\n Number of transactions for training ", size(trans_train, 1), ".")
                print("\n Number of transactions for validation ", size(trans_test, 1), ".")

                # SKU weight
                sku_weight = zeros(Int64, size(trans_train, 2)) .= 1

                #  Start the heuristics and optmisations

                ## Start QMK heuristic with Gurobi as solver
                if start[1, :QMK] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation = MQKP(trans_train, capacity, sku_weight, abort, "Gurobi", show_opt,
                        cpu_cores, allowed_gap, max_nodes, "QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n    QMK: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / gap: ", round(gap_optimisation, digits=6),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="QMK",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "QMK failed" exception=(e, bt)
                        log_failure(logio, "QMK", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start QMK heuristic with Juniper as solver
                if start[1, :QMKJ] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation = MQKP(trans_train, capacity, sku_weight, abort, "Juniper", show_opt,
                        cpu_cores, allowed_gap, max_nodes, "QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n   QMKJ: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / gap: ", round(gap_optimisation, digits=6),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="QMKJ",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "QMKJ failed" exception=(e, bt)
                        log_failure(logio, "QMKJ", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start QMK heuristic with SCIP as solver
                if start[1, :QMKS] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation = MQKP(trans_train, capacity, sku_weight, abort, "scip", show_opt,
                        cpu_cores, allowed_gap, max_nodes, "QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n   QMKS: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / gap: ", round(gap_optimisation, digits=6),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="QMKS",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "QMKS failed" exception=(e, bt)
                        log_failure(logio, "QMKS", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start chi square heuristic without local search
                if start[1, :CHIM] == 1
                    try
                    for sig in sig_levels
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, ls = CHISQUAREHEUR(trans_train, capacity, sig, 0, sku_weight, chistatus)
                        parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                        parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print("\n   CHIM: parcels test data: ", parcels_benchmark,
                            " / parcels training data: ", parcels_train,
                            " / flex: ", flex_benchmark,
                            " / time: ", round(time_benchmark, digits=3),
                            " / local search: ", ls,
                            " / sig: ", sig,
                            " / warehouse: ", sum(W, dims=1))

                        push!(benchmark, (dependency=dependency,
                            skus=skus_benchmark[a],
                            wareh=length(capacity),
                            diff=diff_benchmark[a],
                            buffer=buff_benchmark[a],
                            mode="CHIM_$sig",
                            benchiter=benchnr,
                            orders=size(trans, 1),
                            train_test=train_test,
                            parcel_train=parcels_train,
                            parcel_test=parcels_benchmark,
                            flexibility=flex_benchmark,
                            duration=time_benchmark,
                            cap_used=sum(W),
                            local_search=ls,
                            gap=0))
                    end
                    catch e
                        bt = catch_backtrace()
                        @warn "CHIM failed" exception=(e, bt)
                        log_failure(logio, "CHIM", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start chi square heuristic with local search
                if start[1, :CHI] == 1
                    try
                    for sig in sig_levels
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W, ls = CHISQUAREHEUR(trans_train, capacity, sig, max_ls, sku_weight, chistatus)
                        parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                        parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                        flex_benchmark = FLEXIBILITY(trans_test, W)
                        print("\n    CHI: parcels test data: ", parcels_benchmark,
                            " / parcels training data: ", parcels_train,
                            " / flex: ", flex_benchmark,
                            " / time: ", round(time_benchmark, digits=3),
                            " / local search: ", ls,
                            " / sig: ", sig,
                            " / warehouse: ", sum(W, dims=1))

                        push!(benchmark, (dependency=dependency,
                            skus=skus_benchmark[a],
                            wareh=length(capacity),
                            diff=diff_benchmark[a],
                            buffer=buff_benchmark[a],
                            mode="CHI_$sig",
                            benchiter=benchnr,
                            orders=size(trans, 1),
                            train_test=train_test,
                            parcel_train=parcels_train,
                            parcel_test=parcels_benchmark,
                            flexibility=flex_benchmark,
                            duration=time_benchmark,
                            cap_used=sum(W),
                            local_search=ls,
                            gap=0))
                    end
                    catch e
                        bt = catch_backtrace()
                        @warn "CHI failed" exception=(e, bt)
                        log_failure(logio, "CHI", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start our reproduction of the  K-LINKS heuristic by
                ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
                if start[1, :KL] == 1 && sum(capacity) == length(sku_weight)
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, ls = KLINKS(trans_train, capacity, trials, stagnant, strategy, abort, klinkstatus)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n     KL: parcels train data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / local search: ", ls,
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="KL",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=ls,
                        gap=0))
                    catch e
                        bt = catch_backtrace()
                        @warn "KL failed" exception=(e, bt)
                        log_failure(logio, "KL", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start our reproduction of the  K-LINKS optimization with SBB by
                ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
                if start[1, :KLQ] == 1 && sum(capacity) == length(sku_weight)
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation = MQKP(trans_train, capacity, sku_weight, abort, "Gurobi", show_opt,
                        cpu_cores, allowed_gap, max_nodes, "QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n    KLQ: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / gap: ", round(gap_optimisation, digits=6),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="KLQ",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "KLQ failed" exception=(e, bt)
                        log_failure(logio, "KLQ", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start our reproduction of the greedy orders heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :GO] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = GREEDYORDERS(trans_train, capacity, sku_weight)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n     GO: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="GO",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=0))
                    catch e
                        bt = catch_backtrace()
                        @warn "GO failed" exception=(e, bt)
                        log_failure(logio, "GO", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start our reproduction of the greedy pairs heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :GP] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = GREEDYPAIRS(trans_train, capacity, sku_weight)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n     GP: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="GP",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=0))
                    catch e
                        bt = catch_backtrace()
                        @warn "GP failed" exception=(e, bt)
                        log_failure(logio, "GP", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start our reproduction of the greedy seeds heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :GS] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = GREEDYSEEDS(trans_train, capacity, sku_weight)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n     GS: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="GS",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=0))
                    catch e
                        bt = catch_backtrace()
                        @warn "GS failed" exception=(e, bt)
                        log_failure(logio, "GS", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start our reproduction of the  bestselling heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1, :BS] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = BESTSELLING(trans_train, capacity, sku_weight)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n     BS: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="BS",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=0))
                    catch e
                        bt = catch_backtrace()
                        @warn "BS failed" exception=(e, bt)
                        log_failure(logio, "BS", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start the extended MCI heuristic for D warehouses by
                ## Lin et al. (2025) https://doi.org/10.1111/poms.14114
                if start[1, :EMCI] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = EMCIALLOC(trans_train, capacity, sku_weight)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n   EMCI: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="EMCI",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=0))
                    catch e
                        bt = catch_backtrace()
                        @warn "EMCI failed" exception=(e, bt)
                        log_failure(logio, "EMCI", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start the iterative improvement heuristic with Gurobi by
                ## Lin et al. (2025) https://doi.org/10.1111/poms.14114
                if start[1, :IIH] == 1 && length(capacity) == 2 && sum(capacity) > sum(sku_weight)
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation, ls = IIH(trans_train, capacity, sku_weight,
                        abort, "Gurobi", show_opt, cpu_cores, allowed_gap, max_nodes,
                        max_iih_iterations, epsilon_iih)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n    IIH: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / gap: ", round(gap_optimisation, digits=6),
                        " / iterations: ", ls,
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="IIH",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=ls,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "IIH failed" exception=(e, bt)
                        log_failure(logio, "IIH", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start the iterative improvement heuristic with SCIP by
                ## Lin et al. (2025) https://doi.org/10.1111/poms.14114
                if start[1, :IIHS] == 1 && length(capacity) == 2 && sum(capacity) > sum(sku_weight)
                    try
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation, ls = IIH(trans_train, capacity, sku_weight,
                        abort, "scip", show_opt, cpu_cores, allowed_gap, max_nodes,
                        max_iih_iterations, epsilon_iih)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n   IIHS: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / gap: ", round(gap_optimisation, digits=6),
                        " / iterations: ", ls,
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="IIHS",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=ls,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "IIHS failed" exception=(e, bt)
                        log_failure(logio, "IIHS", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Start the search for optimal solution with the solver CPLEX
                ## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
                ## FULLOPTUEQ if SKUs can be allocated multiple times
                if start[1, :OPT] == 1
                    try
                    sleep(0.01)
                    GC.gc()
                    if sum(capacity) == size(trans, 2)
                        time_benchmark = @elapsed W, gap_optimisation, popt = FULLOPTEQ(trans_train, capacity, abort, show_opt,
                            cpu_cores, allowed_gap, max_nodes)
                    else
                        time_benchmark = @elapsed W, gap_optimisation, popt = FULLOPTUEQ(trans_train, capacity, abort, show_opt,
                            cpu_cores, allowed_gap, max_nodes)
                    end
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    flex_benchmark = FLEXIBILITY(trans_test, W)
                    print("\n    OPT: parcels test data: ", parcels_benchmark,
                        " / parcels training data: ", parcels_train,
                        " / flex: ", flex_benchmark,
                        " / time: ", round(time_benchmark, digits=3),
                        " / warehouse: ", sum(W, dims=1))

                    push!(benchmark, (dependency=dependency,
                        skus=skus_benchmark[a],
                        wareh=length(capacity),
                        diff=diff_benchmark[a],
                        buffer=buff_benchmark[a],
                        mode="OPT",
                        benchiter=benchnr,
                        orders=size(trans, 1),
                        train_test=train_test,
                        parcel_train=parcels_train,
                        parcel_test=parcels_benchmark,
                        flexibility=flex_benchmark,
                        duration=time_benchmark,
                        cap_used=sum(W),
                        local_search=0,
                        gap=gap_optimisation))
                    catch e
                        bt = catch_backtrace()
                        @warn "OPT failed" exception=(e, bt)
                        log_failure(logio, "OPT", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                    end
                end

                ## Benchmark the random allocation of SKUs
                try
                sleep(0.01)
                GC.gc()
                time_benchmark = @elapsed parcels_benchmark, flex_benchmark = RANDOMBENCH(trans_test, capacity, iterations, sku_weight, combination)
                parcels_train, _ = RANDOMBENCH(trans_train, capacity, iterations, sku_weight, combination)
                print("\n    RND: parcels test data: ", parcels_benchmark,
                    " / parcels training data: ", parcels_train,
                    " / flex: ", flex_benchmark,
                    " / time: ", round(time_benchmark, digits=3),
                    " / benchmarks: ", iterations,
                    " / warehouse: ", sum(W, dims=1), "\n")

                push!(benchmark, (dependency=dependency,
                    skus=skus_benchmark[a],
                    wareh=length(capacity),
                    diff=diff_benchmark[a],
                    buffer=buff_benchmark[a],
                    mode="RND",
                    benchiter=benchnr,
                    orders=size(trans, 1),
                    train_test=train_test,
                    parcel_train=parcels_train,
                    parcel_test=parcels_benchmark,
                    flexibility=flex_benchmark,
                    duration=time_benchmark,
                    cap_used=sum(W),
                    local_search=0,
                    gap=0))
                catch e
                    bt = catch_backtrace()
                    @warn "RND failed" exception=(e, bt)
                    log_failure(logio, "RND", e, bt, benchnr, skus_benchmark[a], capacity, orders)
                end

                # Free training/test data before next iteration
                trans_train = nothing
                trans_test = nothing
                sku_weight = nothing
                W = nothing
                GC.gc()

                # Export the results after each stage
                CSV.write("results/$(experiment)_benchmark_$dependency.csv", benchmark)
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
