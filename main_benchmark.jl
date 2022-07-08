function BENCHMARK(capacity_benchmark::Array{Int64,2},
                   skus_benchmark::Vector{Int64},
                   start::DataFrame,
                   orders::Int64,
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
                   sig::Float64)
    
    # Start the benchmark of the data set
    ## Create dataframes for the export of the results
    parcels_benchmark = Array{Float64,2}(undef,size(capacity_benchmark,1), 3+size(start,2))  .= 0

    ### Rename the columns appropriately depending on the selected columns in "start"
    parcels_benchmark = DataFrame(parcels_benchmark, [:wareh, 
                                                      :capacity,
                                                      :buffer, 
                                                      :QMKO,
                                                      :QMK, 
                                                      :CHIM, 
                                                      :CHI, 
                                                      :KL,
                                                      :KLQ,
                                                      :GO,
                                                      :GP,
                                                      :GS,
                                                      :BS,
                                                      :OPT, 
                                                      :RND])
    names_out = String[]
    start[1,:RND] = 1
    for i = 1:size(start,2)
        if start[1,i] == 0
            names_out = vcat(names_out,names(start)[i])
        end
    end
    parcels_benchmark  = select!(parcels_benchmark, Not(names_out))

    ### Copy the correctly named dataframe for the export dataframes
    time_benchmark    = copy(parcels_benchmark)
    parcels_train     = copy(parcels_benchmark)
    gap_optimisation  = copy(parcels_benchmark)

    select!(gap_optimisation, Not(:RND))

    ### Preallocate a transactional data set container
    trans = spzeros(Bool,0,0)

    # Iterate all capacity constellations
    for a = 1:size(capacity_benchmark,1)        
        ## Load the capacity of each individual run
        ### Note: It has to be sorted starting with the largest capacity
        capacity = Array{Int64,1}(undef,count(x -> x > 0, capacity_benchmark[a,:]))
        for b = 1:size(capacity,1)
            capacity[b] = capacity_benchmark[a,b]
        end
        print("\ncapacity constellation: ",a," of ",size(capacity_benchmark,1),
            "\ncapacity: ",capacity)
        
        ## Create all possible capacity combinations for the parcel Benchmark
        combination = COMBINEWAREHOUSES(capacity)

        ## Generate artificial random transactions without dependencies if there is no transactional dataset
        if isfile("transactions/transactions_$experiment.csv") == false
            if  size(trans,2) == skus_benchmark[a]
                    print("\nReused transactions from previous run.")
            else
                print("\nstarting generation of transactions.")
                time = @elapsed trans = RANDOMTRANS(skus_benchmark[a],orders,skus_in_order,sku_frequency,
                                                    ceil(Int64,max(skus_benchmark[a]/50,10)),
                                                    min_dependence,max_dependence,
                                                    group_link,ind_chance,one_direction,
                                                    multi_relatio)
                print("\ntransactions generated after ", round(time,digits = 3)," seconds.")
            end
        end

        ##  Note the numer of warehouses as well as the capacity in the export dataframes
        parcels_benchmark[a,:wareh] = 
        time_benchmark[a,:wareh] = 
        parcels_train[a,:wareh] =
        gap_optimisation[a,:wareh] = size(capacity,1)
        parcels_benchmark[a,:capacity] = 
        time_benchmark[a,:capacity] = 
        parcels_train[a,:capacity] = 
        gap_optimisation[a,:capacity] = sum(capacity)
        parcels_benchmark[a,:buffer] = 
        time_benchmark[a,:buffer] = 
        parcels_train[a,:buffer] =
        gap_optimisation[a,:buffer] = sum(capacity) - size(trans,2)

        #  Split the data into training and test data
        if train_test > 0.00
            cut = round(Int64,size(trans,1) * train_test)
            trans_train = trans[1:cut,:]
            trans_test  = trans[(cut+1):size(trans,1),:]
        else
            trans_train = trans_test = trans
        end
        print("\nNumber of transactions for training ", size(trans_train,1),".") 
        print("\nNumber of transactions for validation ",size(trans_test,1),".")

        #  Start the heuristics and optmisations
        ## Start QMK heuristic to find the optimal solution with the solver CPLEX
        if start[1,:QMKO] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:QMKO] += @elapsed W,gap_optimisation[a,:QMKO] = MQKP(trans_train,capacity,abort,"CPLEX",show_opt,
                                                                                       cpu_cores,allowed_gap,max_nodes,"QMK")
            parcels_benchmark[a,:QMKO] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:QMKO] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n   QMKO: parcels test data: ", parcels_benchmark[a,:QMKO], 
                  " / parcels training data: ", parcels_train[a,:QMKO], 
                  " / time: ",round(time_benchmark[a,:QMKO],digits = 3))
            sleep(0.1)
        end

        ## Start QMK heuristic with SBB as solver
        if start[1,:QMK] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:QMK] += @elapsed W,gap_optimisation[a,:QMK] = MQKP(trans_train,capacity,abort, "SBB",show_opt,
                                                                                 cpu_cores,allowed_gap,max_nodes,"QMK")
            parcels_benchmark[a,:QMK] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:QMK] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n    QMK: parcels test data: ", parcels_benchmark[a,:QMK], 
                  " / parcels training data: ", parcels_train[a,:QMK],  
                  " / time: ", round(time_benchmark[a,:QMK], digits = 3))
            sleep(0.1)
        end

        ## Start chi square heuristic without local search
        if start[1,:CHIM] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:CHIM] += @elapsed W = CHISQUAREHEUR(trans_train,capacity,sig,false,show_opt)
            parcels_benchmark[a,:CHIM] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:CHIM] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n   CHIM: parcels test data: ", parcels_benchmark[a,:CHIM], 
                  " / parcels training data: ", parcels_train[a,:CHIM],  
                  " / time: ",round(time_benchmark[a,:CHIM], digits = 3))
            sleep(0.1)
        end

        ## Start chi square heuristic with local search
        if start[1,:CHI] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:CHI] += @elapsed W = CHISQUAREHEUR(trans_train,capacity,sig,true,show_opt)
            parcels_benchmark[a,:CHI] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:CHI] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n    CHI: parcels test data: ", parcels_benchmark[a,:CHI], 
                  " / parcels training data: ", parcels_train[a,:CHI],  
                  " / time: ",round(time_benchmark[a,:CHI], digits = 3))
            sleep(0.1)
        end

        ## Start our reproduction of the  K-LINKS heuristic by
        ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
        if  start[1,:KL] == 1 #&& sum(capacity) == size(trans,2)
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:KL] += @elapsed W = KLINKS(trans_train,capacity,trials,stagnant,strategy,klinkstatus)
            parcels_benchmark[a,:KL] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:KL] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n     KL: parcels train data: ", parcels_benchmark[a,:KL], 
                  " / parcels training data: ", parcels_train[a,:KL],  
                  " / time: ",round(time_benchmark[a,:KL], digits = 3))
            sleep(0.1)
        end

        ## Start our reproduction of the  K-LINKS optimization with SBB by
        ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
        if  start[1,:KLQ] == 1 #&& sum(capacity) == size(trans,2)
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:KLQ] += @elapsed W, gap_optimisation[a,:KLQ] = MQKP(trans_train,capacity,abort,"SBB",show_opt,
                                                                                        cpu_cores,allowed_gap,max_nodes,"QMK")
            parcels_benchmark[a,:KLQ] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:KLQ] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n    KLQ: parcels test data: ", parcels_benchmark[a,:KLQ], 
                  " / parcels training data: ", parcels_train[a,:KLQ],  
                  " / time: ",round(time_benchmark[a,:KLQ], digits = 3))
            sleep(0.1)
        end

        ## Start our reproduction of the greedy orders heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GO] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:GO] += @elapsed W = GREEDYORDERS(trans_train,capacity)
            parcels_benchmark[a,:GO] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:GO] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n     GO: parcels test data: ", parcels_benchmark[a,:GO], 
                  " / parcels training data: ", parcels_train[a,:GO],  
                  " / time: ",round(time_benchmark[a,:GO], digits = 3))
            sleep(0.1)
        end

        ## Start our reproduction of the greedy pairs heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GP] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:GP] += @elapsed W = GREEDYPAIRS(trans_train,capacity)
            parcels_benchmark[a,:GP] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:GP] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n     GP: parcels test data: ", parcels_benchmark[a,:GP], 
                  " / parcels training data: ", parcels_train[a,:GP],  
                  " / time: ",round(time_benchmark[a,:GP], digits = 3))
            sleep(0.1)
        end

        ## Start our reproduction of the greedy seeds heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GS] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:GS] += @elapsed W = GREEDYSEEDS(trans_train,capacity)
            parcels_benchmark[a,:GS] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:GS] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n     GS: parcels test data: ", parcels_benchmark[a,:GS], 
                  " / parcels training data: ", parcels_train[a,:GS],  
                  " / time: ",round(time_benchmark[a,:GS], digits = 3))
            sleep(0.1)
        end

        ## Start our reproduction of the  bestselling heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if  start[1,:BS] == 1
            sleep(0.1)
            GC.gc()
            time_benchmark[a,:BS] += @elapsed W = BESTSELLING(trans_train,capacity)
            parcels_benchmark[a,:BS] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:BS] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n     BS: parcels test data: ", parcels_benchmark[a,:BS], 
                  " / parcels training data: ", parcels_train[a,:BS],  
                  " / time: ",round(time_benchmark[a,:BS], digits = 3))
            sleep(0.1)
        end

        ## Start the search for optimal solution with the solver CPLEX
        ## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
        ## FULLOPTUEQ if SKUs can be allocated multiple times
        if start[1,:OPT] == 1
            sleep(0.1)
            GC.gc()
            if sum(capacity) == size(trans,2)
                time_benchmark[a,:OPT] += @elapsed W,gap_optimisation[a,:OPT],popt = FULLOPTEQ(trans_train,capacity,abort,show_opt,
                                                                                            cpu_cores,allowed_gap,max_nodes)
            else
                time_benchmark[a,:OPT] += @elapsed W,gap_optimisation[a,:OPT],popt = FULLOPTUEQ(trans_train,capacity,abort,show_opt,
                                                                                             cpu_cores,allowed_gap,max_nodes)
            end
            parcels_benchmark[a,:OPT] = PARCELSSEND(trans_test, W, capacity, combination)
            parcels_train[a,:OPT] = PARCELSSEND(trans_train, W, capacity, combination)
            print("\n    OPT: parcels test data: ", parcels_benchmark[a,:OPT], 
                  " / parcels training data: ", parcels_train[a,:OPT],  
                  " / time: ",round(time_benchmark[a,:OPT], digits = 3))
            sleep(0.1)
        end

        ## Benchmark the random allocation of SKUs
        sleep(0.1)
        time_benchmark[a,:RND] += @elapsed parcels_benchmark[a,:RND] = RANDOMBENCH(trans_test,capacity,iterations,combination)
        parcels_train[a,:RND] = RANDOMBENCH(trans_train,capacity,iterations,combination)
        print("\n    RND: parcels test data: ", parcels_benchmark[a,:RND],
              " / parcels training data: ", parcels_train[a,:RND],"\n")
        sleep(0.1)

        # Export the results after each stage
        CSV.write("results/$(experiment)_a_parcels_sent_$dependency.csv",       parcels_benchmark)
        CSV.write("results/$(experiment)_b_duration_$dependency.csv",           time_benchmark)
        CSV.write("results/$(experiment)_c_parcels_train_$dependency.csv",      parcels_train)
        CSV.write("results/$(experiment)_d_optimisation_gap_$dependency.csv",   gap_optimisation)
    end
    print("\n### Final Report ###",
          "\nparcels send: \n",parcels_benchmark,
          "\nTime needed: \n",time_benchmark, "\n\n")
    
    return parcels_benchmark::DataFrame, 
           time_benchmark::DataFrame, 
           parcels_train::DataFrame, 
           gap_optimisation::DataFrame
end