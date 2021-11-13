function BENCHMARK(capacity_benchmark::Array{Int64,2},
                   start::DataFrame,
                   orders::Int64,
                   trials::Int64,
                   stagnant::Int64,
                   strategy::Int64,
                   klinkstatus::Int64,
                   abort::Int64,
                   iterations::Int64,
                   show_opt::Int64,
                   cpu_cores::Int64,
                   allowed_gap::Float64,
                   max_nodes::Int64,
                   sig::Float64)
    
    # Start the benchmark of the data set
    ## Create dataframes for the export of the results
    parcels_benchmark = Array{Float64,2}(undef,size(capacity_benchmark,1), 2+size(start,2))  .= 0

    ### Rename the columns appropriately depending on the selected columns in "start"
    parcels_benchmark = DataFrame(parcels_benchmark, [:wareh, 
                                                      :capacity, 
                                                      :QMKOPT, 
                                                      :QMK, 
                                                      :CHISOL, 
                                                      :CHI, 
                                                      :KLINK, 
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
    cap_used          = copy(parcels_benchmark)
    parcel_reduction  = copy(parcels_benchmark)
    split_reduction   = copy(parcels_benchmark)
    gap_optimisation  = copy(parcels_benchmark)

    ## import the transactional data and display the strength of the correlation within the transactional dataset
    if isfile("transactions/transactions_$experiment.csv") == true
        trans = readdlm("transactions/transactions_$experiment.csv", ',', Int64)
        corrtrans = cor(trans, dims=1)
        display(heatmap(corrtrans, 
                    yaxis=("SKU"), 
                    xaxis=("SKU"), 
                    dpi = 300, 
                    c = cgrad(:bone_1), 
                    clim = (-1,1)))
                    savefig("graphs/cor_$experiment.pdf")
    end

    # Iterate all capacity constellations
    print("\n\n### Benchmark started ###")
    for a = 1:size(capacity_benchmark,1)        
        ## Load the capacity of each individual run
        ### Note: It has to be sorted starting with the largest capacity
        capacity = Array{Int64,1}(undef,count(x -> x > 0, capacity_benchmark[a,:]))
        for b = 1:size(capacity,1)
            capacity[b] = capacity_benchmark[a,b]
        end
        print("\ncapacity constellation: ",a," of ",size(capacity_benchmark,1),
            "\ncapacity: ",capacity,"\n")
        ##  Note the numer of warehouses as well as the capacity in the export dataframes
        parcels_benchmark[a,:wareh] = 
        time_benchmark[a,:wareh] = 
        cap_used[a,:wareh] =
        parcel_reduction[a,:wareh] = 
        split_reduction[a,:wareh] = 
        gap_optimisation[a,:wareh] = size(capacity,1)
        parcels_benchmark[a,:capacity] = 
        time_benchmark[a,:capacity] = 
        cap_used[a,:capacity] = 
        parcel_reduction[a,:capacity] = 
        split_reduction[a,:capacity] =
        gap_optimisation[a,:capacity] = sum(capacity)

        ## Create all possible capacity combinations for the parcel Benchmark
        combination = COMBINEWAREHOUSES(capacity)

        ## Generate artificial random transactions without dependencies if there is no transactional dataset
        if isfile("transactions/transactions_$experiment.csv") == false
            trans = RANDOMIND(orders,sum(capacity_benchmark[a,:]))
        end

        #  Start the heuristics and optmisations
        ## Start QMK heuristic to find the optimal solution with the solver CPLEX
        if start[1,:QMKOPT] == 1
            time_benchmark[a,:QMKOPT]  = @elapsed Q = COAPPEARENCE(trans::Array{Int64,2})
            time_benchmark[a,:QMKOPT] += @elapsed W,gap_optimisation[a,:QMKOPT] = MQKP(trans::Array{Int64,2},
                                                                                       capacity::Array{Int64,1},
                                                                                       Q::Array{Int64,2},
                                                                                       abort::Int64,
                                                                                       "CPLEX",
                                                                                       show_opt::Int64,
                                                                                       cpu_cores::Int64,
                                                                                       allowed_gap::Float64,
                                                                                       max_nodes::Int64)
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:QMKOPT] = sum(W)
            parcels_benchmark[a,:QMKOPT] = sum(parcel)
            print("\n      qmkopt: parcels after optimisation: ", parcels_benchmark[a,:QMKOPT], 
                  " / capacity_used: ", cap_used[a,:QMKOPT], " / time: ",time_benchmark[a,:QMKOPT])
        end
        
        ## Start QMK heuristic with BONMIN as solver
        if start[1,:QMK] == 1
            time_benchmark[a,:QMK]  = @elapsed Q = COAPPEARENCE(trans::Array{Int64,2})
            time_benchmark[a,:QMK] += @elapsed W,gap_optimisation[a,:QMK] = MQKP(trans::Array{Int64,2},
                                                                                 capacity::Array{Int64,1},
                                                                                 Q::Array{Int64,2},
                                                                                 abort::Int64,
                                                                                 "BONMIN",
                                                                                 show_opt::Int64,
                                                                                 cpu_cores::Int64,
                                                                                 allowed_gap::Float64,
                                                                                 max_nodes::Int64)
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:QMK] = sum(W)
            parcels_benchmark[a,:QMK] = sum(parcel)
            print("\n         mqk: parcels after optimisation: ", parcels_benchmark[a,:QMK], 
                  " / capacity_used: ", cap_used[a,:QMK],  " / time: ",time_benchmark[a,:QMK])
        end

        ## Start CHI square heuristic without local search
        if start[1,:CHISOL] == 1
            time_benchmark[a,:CHISOL]  = @elapsed Q = COAPPEARENCE(trans::Array{Int64,2})
            time_benchmark[a,:CHISOL] += @elapsed W = CHISQUAREHEUR(trans::Array{Int64,2},
                                                                    capacity::Array{Int64,1},
                                                                    Q::Array{Int64,2},
                                                                    sig::Float64) 
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:CHISOL] = sum(W)
            parcels_benchmark[a,:CHISOL] = sum(parcel)
            print("\n      chisol: parcels after optimisation: ", parcels_benchmark[a,:CHISOL], 
                  " / capacity_used: ", cap_used[a,:CHISOL],  " / time: ",time_benchmark[a,:CHISOL])
        end

        ## Start CHI square heuristic with local search
        if start[1,:CHI] == 1
            time_benchmark[a,:CHI]  = @elapsed Q = COAPPEARENCE(trans::Array{Int64,2})
            time_benchmark[a,:CHI] += @elapsed W = CHISQUAREHEUR(trans::Array{Int64,2},
                                                                 capacity::Array{Int64,1},
                                                                 Q::Array{Int64,2},
                                                                 sig::Float64)
            time_benchmark[a,:CHI] += @elapsed W = LOCALSEARCH(W::Array{Int64,2}, Q::Array{Int64,2})
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:CHI] = sum(W)
            parcels_benchmark[a,:CHI] = sum(parcel)
            print("\n         chi: parcels after optimisation: ", parcels_benchmark[a,:CHI], 
                  " / capacity_used: ", cap_used[a,:CHI],  " / time: ",time_benchmark[a,:CHI])
        end

        ## Start our reproduction of the  K-LINKS heuristic by
        ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
        if  start[1,:KLINK] == 1 && sum(capacity) == size(trans,2)
            time_benchmark[a,:KLINK] = @elapsed W = KLINKS(trans::Array{Int64,2},
                                                           capacity::Array{Int64,1},
                                                           trials::Int64,
                                                           stagnant::Int64,
                                                           strategy::Int64,
                                                           klinkstatus::Int64)
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:KLINK] = sum(W)
            parcels_benchmark[a,:KLINK] = sum(parcel)
            print("\n     k-links: parcels after optimisation: ", parcels_benchmark[a,:KLINK], 
                  " / capacity_used: ", cap_used[a,:KLINK],  " / time: ",time_benchmark[a,:KLINK])
        end

        ## Start our reproduction of the greedy pairs heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GP] == 1
            time_benchmark[a,:GP] = @elapsed W = GREEDYPAIRS(trans::Array{Int64,2}, 
                                                             capacity::Array{Int64,1})
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:GP] = sum(W)
            parcels_benchmark[a,:GP] = sum(parcel)
            print("\ngreedy pairs: parcels after optimisation: ", parcels_benchmark[a,:GP], 
                  " / capacity_used: ", cap_used[a,:GP],  " / time: ",time_benchmark[a,:GP])
        end

        ## Start our reproduction of the greedy seeds heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GS] == 1
            time_benchmark[a,:GS] = @elapsed W = GREEDYSEEDS(trans::Array{Int64,2},
                                                             capacity::Array{Int64,1})
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:GS] = sum(W)
            parcels_benchmark[a,:GS] = sum(parcel)
            print("\ngreedy seeds: parcels after optimisation: ", parcels_benchmark[a,:GS], 
                  " / capacity_used: ", cap_used[a,:GS],  " / time: ",time_benchmark[a,:GS])
        end

        ## Start our reproduction of the  bestselling heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if  start[1,:BS] == 1  && sum(capacity) > size(trans,2)
            time_benchmark[a,:BS] = @elapsed W = BESTSELLING(trans::Array{Int64,2},
                                                             capacity::Array{Int64,1})
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:BS] = sum(W)
            parcels_benchmark[a,:BS] = sum(parcel)
            print("\n bestselling: parcels after optimisation: ", parcels_benchmark[a,:BS], 
                  " / capacity_used: ", cap_used[a,:BS],  " / time: ",time_benchmark[a,:BS])
        end

        ## Start the search for optimal solution with the solver CPLEX
        ## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
        ## FULLOPTUEQ if SKUs can be allocated multiple times
        if start[1,:OPT] == 1
            if sum(capacity) == size(trans,2)
                time_benchmark[a,:OPT] = @elapsed W,gap_optimisation[a,5],popt = FULLOPTEQ(trans::Array{Int64,2},
                                                                                           capacity::Array{Int64,1},
                                                                                           abort::Int64,
                                                                                           show_opt::Int64,
                                                                                           cpu_cores::Int64,
                                                                                           allowed_gap::Float64,
                                                                                           max_nodes::Int64)
            else
                time_benchmark[a,:OPT] = @elapsed W,gap_optimisation[a,5],popt = FULLOPTUEQ(trans::Array{Int64,2},
                                                                                            capacity::Array{Int64,1},
                                                                                            abort::Int64,
                                                                                            show_opt::Int64,
                                                                                            cpu_cores::Int64,
                                                                                            allowed_gap::Float64,
                                                                                            max_nodes::Int64)
            end
            parcel = PARCELSSEND(trans::Array{Int64,2}, 
                                 W::Array{Int64,2}, 
                                 capacity::Array{Int64,1},
                                 combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:OPT] = sum(W)
            if popt == sum(parcel)
                parcels_benchmark[a,:OPT] = sum(parcel)
                print("\n Solution of optimisation: ",popt, "/ Solution of parcel simulation: ", sum(parcel),
                " / time: ",time_benchmark[a,:OPT])
            else
                print("\n Optimisation Value: ",popt,"/ Parcel Simulation: ",sum(parcel))
                error("Optimisation not working correctly.")
            end
        end

        ## Benchmark the random allocation of SKUs
        time_benchmark[a,:RND] = @elapsed parcels_benchmark[a,:RND] = RANDOMBENCH(trans::Array{Int64,2},
                                                                                  capacity::Array{Int64,1},
                                                                                  iterations::Int64,
                                                                                  combination::Array{Array{Array{Int64,1},1},1})
        print("\n      random: parcels after optimisation: ", parcels_benchmark[a,:RND])

        ## Calculate number of split deliveries
        split_reduction[a:a,3:end] .= parcels_benchmark[a:a,3:end] .- size(trans,1)
        print("\n\n")
    end

    ## Caclculate the improvements of the heuristics compared to the random allocations
    for i = 1:size(parcels_benchmark,1)
        for j = 3:size(parcels_benchmark,2)
            if parcels_benchmark[i,j] > 0
                parcel_reduction[i,j] = (parcels_benchmark[i,:RND]-parcels_benchmark[i,j])/parcels_benchmark[i,:RND]
                split_reduction[i,j] = (split_reduction[i,:RND]-split_reduction[i,j])/(split_reduction[i,:RND])
            else
                parcel_reduction[i,j] = 0
                split_reduction[i,j] = 0
            end
        end
    end

    print("\n### Final Report ###",
          "\nparcels send: \n",parcels_benchmark,
          "\nTime needed: \n",time_benchmark,
          "\nCap used: \n",cap_used,
          "\nParcel Reduction: \n",parcel_reduction,
          "\nSplit Reduction: \n",split_reduction)
    
    return parcels_benchmark::DataFrame, 
           time_benchmark::DataFrame, 
           cap_used::DataFrame, 
           parcel_reduction::DataFrame, 
           split_reduction::DataFrame, 
           gap_optimisation::DataFrame
end