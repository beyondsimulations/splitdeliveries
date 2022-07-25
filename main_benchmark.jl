function BENCHMARK(capacity_benchmark::Array{Int64,2},
                   skus_benchmark::Vector{Int64},
                   start::DataFrame,
                   order_sets::Vector{Int64},
                   max_dependence::Float64,
                   trials::Int64,
                   stagnant::Int64,
                   strategy::Int64,
                   klinkstatus::Int64,
                   abort::Int64,
                   iterations::Int64,
                   show_opt::Bool,
                   cpu_cores::Int64,
                   allowed_gap::Float64,
                   max_nodes::Int64,
                   sig_levels::Vector{Float64},
                   max_ls::Int64,
                   chistatus::Bool,
                   benchiterations::Int64)
    
    # Start the benchmark of the data set
    ## Create dataframes for the export of the results
    benchmark = DataFrame(dependency = String[], 
                          skus = Int64[], 
                          wareh= Int64[], 
                          diff = Float64[], 
                          buffer = Float64[], 
                          mode = String[],
                          benchiter = Int64[],
                          orders = Int64[],
                          train_test = Float64[],
                          parcel_train = Int64[], 
                          parcel_test = Int64[], 
                          duration = Float64[],
                          cap_used = Int64[],
                          local_search = Int64[],
                          gap = Float64[])

    ### Preallocate a transactional data set container
    trans = spzeros(Bool,0,0)

    # Repeat the benchmark benchiter times
    for benchnr = 1:benchiterations
        print("\nStarting benchmark iteration ",benchnr," of ", benchiterations)

        # Iterate all capacity constellations
        for a in axes(capacity_benchmark,1)        
            ## Load the capacity of each individual run
            ### Note: It has to be sorted starting with the largest capacity
            capacity = Array{Int64,1}(undef,count(x -> x > 0, capacity_benchmark[a,:]))
            for b in axes(capacity,1)
                capacity[b] = capacity_benchmark[a,b]
            end
            print("\n Benchmark ", benchnr, " of ", benchiterations,
                 "/ Constellation: ",a," of ",size(capacity_benchmark,1),
                "\n Capacity: ",capacity)
            
            ## Create all possible capacity combinations for the parcel Benchmark
            combination = COMBINEWAREHOUSES(capacity)

            ## Iterate over different number of orders
            for orders in order_sets

                ## Generate artificial random transactions without dependencies if there is no transactional dataset
                if isfile("transactions/transactions_$experiment.csv") == false
                    if  size(trans,2) == skus_benchmark[a] && size(trans,1) == orders
                            print("\n Reused transactions from previous run.")
                    else
                        print("\n Starting generation of transactions.")
                        time = @elapsed trans = RANDOMTRANS(skus_benchmark[a],orders,skus_in_order,sku_frequency,
                                                            ceil(Int64,max(skus_benchmark[a]/50,10)),
                                                            min_dependence,max_dependence,
                                                            group_link,ind_chance,one_direction,
                                                            multi_relatio)
                        print("\n Transactions generated after ", round(time,digits = 3)," seconds.")
                    end
                end

                #  Split the data into training and test data
                if train_test > 0.00
                    cut = round(Int64,size(trans,1) * train_test)
                    trans_train = trans[1:cut,:]
                    trans_test  = trans[(cut+1):size(trans,1),:]
                else
                    trans_train = trans_test = trans
                end
                print("\n Number of transactions for training ", size(trans_train,1),".") 
                print("\n Number of transactions for validation ",size(trans_test,1),".")

                #  Start the heuristics and optmisations
                ## Start QMK heuristic to find the optimal solution with the solver CPLEX
                if start[1,:QMKO] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W,gap_optimisation = MQKP(trans_train,capacity,abort,"CPLEX",show_opt,
                                                                        cpu_cores,allowed_gap,max_nodes,"QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n   QMKO: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train, 
                        " / time: ",round(time_benchmark,digits = 3),
                        " / gap: ",round(gap_optimisation,6))
                    
                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "QMKO",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = gap_optimisation))
                end

                ## Start QMK heuristic with SBB as solver
                if start[1,:QMK] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W,gap_optimisation = MQKP(trans_train,capacity,abort, "SBB",show_opt,
                                                                        cpu_cores,allowed_gap,max_nodes,"QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n    QMK: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ", round(time_benchmark, digits = 3),
                        " / gap: ",round(gap_optimisation,6))

                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "QMK",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = gap_optimisation))
                end

                ## Start chi square heuristic without local search
                if start[1,:CHIM] == 1
                    for sig in sig_levels
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W,ls = CHISQUAREHEUR(trans_train,capacity,sig,0,chistatus)
                        parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                        parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                        print("\n   CHIM: parcels test data: ", parcels_benchmark, 
                            " / parcels training data: ", parcels_train,  
                            " / time: ",round(time_benchmark, digits = 3),
                            " / local search: ", ls,
                            " / sig: ", sig)

                        push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "CHIM_$sig",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = ls,
                                        gap = 0))
                    end 
                end

                ## Start chi square heuristic with local search
                if start[1,:CHI] == 1
                    for sig in sig_levels
                        sleep(0.01)
                        GC.gc()
                        time_benchmark = @elapsed W,ls = CHISQUAREHEUR(trans_train,capacity,sig,max_ls,chistatus)
                        parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                        parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                        print("\n    CHI: parcels test data: ", parcels_benchmark, 
                            " / parcels training data: ", parcels_train,  
                            " / time: ",round(time_benchmark, digits = 3),
                            " / local search: ", ls,
                            " / sig: ", sig)

                        push!(benchmark, (dependency = dependency, 
                                            skus = skus_benchmark[a],  
                                            wareh = length(capacity), 
                                            diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                            buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                            mode = "CHI_$sig",
                                            benchiter = benchnr,
                                            orders = size(trans,1),
                                            train_test = train_test,
                                            parcel_train = parcels_train, 
                                            parcel_test = parcels_benchmark, 
                                            duration = time_benchmark,
                                            cap_used = sum(W),
                                            local_search = ls,
                                            gap = 0))
                    end
                end

                ## Start our reproduction of the  K-LINKS heuristic by
                ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
                if  start[1,:KL] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W,ls = KLINKS(trans_train,capacity,trials,stagnant,strategy,klinkstatus)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n     KL: parcels train data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ", round(time_benchmark, digits = 3),
                        " / local search: ", ls)

                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "KL",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = ls,
                                        gap = 0))
                end

                ## Start our reproduction of the  K-LINKS optimization with SBB by
                ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
                if  start[1,:KLQ] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W, gap_optimisation = MQKP(trans_train,capacity,abort,"SBB",show_opt,
                                                                        cpu_cores,allowed_gap,max_nodes,"QMK")
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n    KLQ: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ", round(time_benchmark, digits = 3),
                        " / gap: ", round(gap_optimisation,6))

                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "KLQ",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = gap_optimisation))
                end

                ## Start our reproduction of the greedy orders heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1,:GO] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = GREEDYORDERS(trans_train,capacity)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n     GO: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))

                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "GO",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = 0))
                end

                ## Start our reproduction of the greedy pairs heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1,:GP] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = GREEDYPAIRS(trans_train,capacity)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n     GP: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))
                    
                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "GP",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = 0))
                end

                ## Start our reproduction of the greedy seeds heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if start[1,:GS] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = GREEDYSEEDS(trans_train,capacity)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n     GS: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))

                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "GS",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = 0))
                end

                ## Start our reproduction of the  bestselling heuristic by
                ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
                if  start[1,:BS] == 1
                    sleep(0.01)
                    GC.gc()
                    time_benchmark = @elapsed W = BESTSELLING(trans_train,capacity)
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n     BS: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))
                        
                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "BS",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = 0))
                end

                ## Start the search for optimal solution with the solver CPLEX
                ## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
                ## FULLOPTUEQ if SKUs can be allocated multiple times
                if start[1,:OPT] == 1
                    sleep(0.01)
                    GC.gc()
                    if sum(capacity) == size(trans,2)
                        time_benchmark = @elapsed W,gap_optimisation,popt = FULLOPTEQ(trans_train,capacity,abort,show_opt,
                                                                                        cpu_cores,allowed_gap,max_nodes)
                    else
                        time_benchmark = @elapsed W,gap_optimisation,popt = FULLOPTUEQ(trans_train,capacity,abort,show_opt,
                                                                                        cpu_cores,allowed_gap,max_nodes)
                    end
                    parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                    parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                    print("\n    OPT: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))
                    
                    push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "OPT",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = gap_optimisation))
                end

                ## Benchmark the random allocation of SKUs
                sleep(0.01)
                time_benchmark = @elapsed parcels_benchmark = RANDOMBENCH(trans_test,capacity,iterations,combination)
                parcels_train = RANDOMBENCH(trans_train,capacity,iterations,combination)
                print("\n    RND: parcels test data: ", parcels_benchmark,
                    " / parcels training data: ", parcels_train,
                    " / time: ", round(time_benchmark, digits = 3),
                    " / benchmarks: ", iterations, "\n")

                push!(benchmark, (dependency = dependency, 
                                        skus = skus_benchmark[a],  
                                        wareh = length(capacity), 
                                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                                        buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                        mode = "RND",
                                        benchiter = benchnr,
                                        orders = size(trans,1),
                                        train_test = train_test,
                                        parcel_train = parcels_train, 
                                        parcel_test = parcels_benchmark, 
                                        duration = time_benchmark,
                                        cap_used = sum(W),
                                        local_search = 0,
                                        gap = 0))

                # Export the results after each stage
                CSV.write("results/$(experiment)_benchmark_$dependency.csv", benchmark)
            end
        end
    end

    print("\n### Finished Benchmark ###")
    
    return benchmark::DataFrame
end