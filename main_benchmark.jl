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
                                                      :QMKLOC, 
                                                      :CHI, 
                                                      :CHILOC, 
                                                      :KLINK, 
                                                      :GP,
                                                      :GPLOC, 
                                                      :GS,
                                                      :GSLOC, 
                                                      :BS,
                                                      :BSLOC, 
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
        trans = trans .* 1.0
        trans = sparse(trans)
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
            "\ncapacity: ",capacity)
            
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
            print("\nstarting generation of transactions.")
            time = @elapsed trans = RANDOMTRANS(skus_benchmark[a],orders,ceil(Int64,skus_benchmark[a]/10),
                                                min_dependence,max_dependence,
                                                group_link,ind_chance,one_direction,
                                                multi_relatio)
            print("\ntransactions generated after ", round(time,digits = 3)," seconds.")
        end

        #  Calculate the coappearance matrix
            time_benchmark[a,3:end] .= coapp_time =  @elapsed Q = COAPPEARENCE(trans)
            if start[1,:KLINK] == 1
                time_benchmark[a,:KLINK] = 0
            end
            if start[1,:OPT] == 1
                time_benchmark[a,:OPT] = 0
            end
            if start[1,:RND] == 1
                time_benchmark[a,:RND] = 0
            end
            print("\ncoappearance matrix: ",round(coapp_time, digits = 3)," seconds\n")

        #  Start the heuristics and optmisations
        ## Start QMK heuristic to find the optimal solution with the solver CPLEX
        if start[1,:QMKOPT] == 1
            time_benchmark[a,:QMKOPT] += @elapsed W,gap_optimisation[a,:QMKOPT] = MQKP(capacity::Array{Int64,1},
                                                                                       Q::Array{Int64,2},
                                                                                       abort::Int64,
                                                                                       "CPLEX",
                                                                                       show_opt::Int64,
                                                                                       cpu_cores::Int64,
                                                                                       allowed_gap::Float64,
                                                                                       max_nodes::Int64)
            parcel = PARCELSSEND(trans, W::Array{Int64,2}, capacity::Array{Int64,1}, combination::Array{Array{Array{Int64,1},1},1})
            cap_used[a,:QMKOPT] = sum(W)
            parcels_benchmark[a,:QMKOPT] = PARCELSSEND(trans, W, capacity, combination)
            print("\n      qmkopt: parcels after optimisation: ", parcels_benchmark[a,:QMKOPT], 
                  " / capacity_used: ", cap_used[a,:QMKOPT], " / time: ",round(time_benchmark[a,:QMKOPT],digits = 3))
        end

        ## Start QMK heuristic with SBB as solver
        if start[1,:QMK] == 1
            time_benchmark[a,:QMK] += @elapsed W,gap_optimisation[a,:QMK] = MQKP(capacity::Array{Int64,1},
                                                                                 Q::Array{Int64,2},
                                                                                 abort::Int64,
                                                                                 "SBB",
                                                                                 show_opt::Int64,
                                                                                 cpu_cores::Int64,
                                                                                 allowed_gap::Float64,
                                                                                 max_nodes::Int64)
            cap_used[a,:QMK] = sum(W)
            parcels_benchmark[a,:QMK] = PARCELSSEND(trans, W, capacity, combination)
            print("\n         mqk: parcels after optimisation: ", parcels_benchmark[a,:QMK], 
                  " / capacity_used: ", cap_used[a,:QMK],  " / time: ", round(time_benchmark[a,:QMK], digits = 3))
        end

        ## Start QMK heuristic with BONMIN as solver
        if start[1,:QMKLOC] == 1
            if start[1,:QMK] == 1
                time_benchmark[a,:QMK] = time_benchmark[a,:QMKLOC]
            else
                time_benchmark[a,:QMKLOC] += @elapsed W,gap_optimisation[a,:QMKLOC] = MQKP(capacity::Array{Int64,1},
                                                                                            Q::Array{Int64,2},
                                                                                            abort::Int64,
                                                                                            "SBB",
                                                                                            show_opt::Int64,
                                                                                            cpu_cores::Int64,
                                                                                            allowed_gap::Float64,
                                                                                            max_nodes::Int64)
            end
            time_benchmark[a,:QMKLOC] += @elapsed W = LOCALSEARCH(W::Array{Int64,2}, Q::Array{Int64,2})
            cap_used[a,:QMKLOC] = sum(W)
            parcels_benchmark[a,:QMKLOC] = PARCELSSEND(trans, W, capacity, combination)
            print("\n     mqk+loc: parcels after optimisation: ", parcels_benchmark[a,:QMKLOC], 
                  " / capacity_used: ", cap_used[a,:QMKLOC],  " / time: ", round(time_benchmark[a,:QMKLOC], digits = 3))
        end

        ## Start CHI square heuristic without local search
        if start[1,:CHI] == 1
            time_benchmark[a,:CHI] += @elapsed W = CHISQUAREHEUR(trans, capacity::Array{Int64,1},
                                                                 Q::Array{Int64,2},
                                                                 sig::Float64)
            cap_used[a,:CHI] = sum(W)
            parcels_benchmark[a,:CHI] = PARCELSSEND(trans, W, capacity, combination)
            print("\n         chi: parcels after optimisation: ", parcels_benchmark[a,:CHI], 
                  " / capacity_used: ", cap_used[a,:CHI],  " / time: ",round(time_benchmark[a,:CHI], digits = 3))
        end

        ## Start CHI square heuristic with local search
        if start[1,:CHILOC] == 1
            if start[1,:CHI] == 1
                time_benchmark[a,:CHILOC] = time_benchmark[a,:CHI]
            else
                time_benchmark[a,:CHILOC] += @elapsed W = CHISQUAREHEUR(trans, capacity::Array{Int64,1},
                                                                        Q::Array{Int64,2},
                                                                        sig::Float64)
            end
            time_benchmark[a,:CHILOC] += @elapsed W = LOCALSEARCHCHI(W::Array{Int64,2}, Q::Array{Int64,2})
            cap_used[a,:CHILOC] = sum(W)
            parcels_benchmark[a,:CHILOC] = PARCELSSEND(trans, W, capacity, combination)
            print("\n     chi+loc: parcels after optimisation: ", parcels_benchmark[a,:CHILOC], 
                  " / capacity_used: ", cap_used[a,:CHILOC],  " / time: ",round(time_benchmark[a,:CHILOC], digits = 3))
        end

        ## Start our reproduction of the  K-LINKS heuristic by
        ## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
        if  start[1,:KLINK] == 1 && sum(capacity) == size(trans,2)
            time_benchmark[a,:KLINK] += @elapsed W = KLINKS(trans::SparseMatrixCSC{Float64, Int64},
                                                           capacity::Array{Int64,1},
                                                           trials::Int64,
                                                           stagnant::Int64,
                                                           strategy::Int64,
                                                           klinkstatus::Int64)
            cap_used[a,:KLINK] = sum(W)
            parcels_benchmark[a,:KLINK] = PARCELSSEND(trans, W, capacity, combination)
            print("\n     k-links: parcels after optimisation: ", parcels_benchmark[a,:KLINK], 
                  " / capacity_used: ", cap_used[a,:KLINK],  " / time: ",round(time_benchmark[a,:KLINK], digits = 3))
        end

        ## Start our reproduction of the greedy pairs heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GP] == 1
            time_benchmark[a,:GP] += @elapsed W = GREEDYPAIRS(Q::Matrix{Int64},
                                                              capacity::Array{Int64,1})
            cap_used[a,:GP] = sum(W)
            parcels_benchmark[a,:GP] = PARCELSSEND(trans, W, capacity, combination)
            print("\n          gp: parcels after optimisation: ", parcels_benchmark[a,:GP], 
                  " / capacity_used: ", cap_used[a,:GP],  " / time: ",round(time_benchmark[a,:GP], digits = 3))
        end

        ## Start our reproduction of the greedy pairs heuristic with an additional local search by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GPLOC] == 1
            time_benchmark[a,:GPLOC] += @elapsed W = GREEDYPAIRS(Q::Matrix{Int64},
                                                                 capacity::Array{Int64,1})
            time_benchmark[a,:GPLOC] += @elapsed W = LOCALSEARCH(W::Array{Int64,2}, Q::Array{Int64,2})
            cap_used[a,:GPLOC] = sum(W)
            parcels_benchmark[a,:GPLOC] = PARCELSSEND(trans, W, capacity, combination)
            print("\n      gp+loc: parcels after optimisation: ", parcels_benchmark[a,:GPLOC], 
                  " / capacity_used: ", cap_used[a,:GPLOC],  " / time: ",round(time_benchmark[a,:GPLOC], digits = 3))
        end

        ## Start our reproduction of the greedy seeds heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GS] == 1
            time_benchmark[a,:GS] += @elapsed W = GREEDYSEEDS(trans, Q::Matrix{Int64},
                                                             capacity::Array{Int64,1})
            cap_used[a,:GS] = sum(W)
            parcels_benchmark[a,:GS] = PARCELSSEND(trans, W, capacity, combination)
            print("\n          gs: parcels after optimisation: ", parcels_benchmark[a,:GS], 
                  " / capacity_used: ", cap_used[a,:GS],  " / time: ",round(time_benchmark[a,:GS], digits = 3))
        end

        ## Start our reproduction of the greedy seeds heuristic with an additional local search by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if start[1,:GSLOC] == 1
            time_benchmark[a,:GSLOC] += @elapsed W = GREEDYSEEDS(trans, Q::Matrix{Int64},
                                                                 capacity::Array{Int64,1})
            time_benchmark[a,:GSLOC] += @elapsed W = LOCALSEARCH(W::Array{Int64,2}, Q::Array{Int64,2})
            cap_used[a,:GSLOC] = sum(W)
            parcels_benchmark[a,:GSLOC] = PARCELSSEND(trans, W, capacity, combination)
            print("\n      gs+loc: parcels after optimisation: ", parcels_benchmark[a,:GSLOC], 
                  " / capacity_used: ", cap_used[a,:GSLOC],  " / time: ",round(time_benchmark[a,:GSLOC], digits = 3))
        end

        ## Start our reproduction of the  bestselling heuristic by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if  start[1,:BS] == 1
            time_benchmark[a,:BS] += @elapsed W = BESTSELLING(trans, Q::Matrix{Int64},
                                                              capacity::Array{Int64,1})
            cap_used[a,:BS] = sum(W)
            parcels_benchmark[a,:BS] = PARCELSSEND(trans, W, capacity, combination)
            print("\n          bs: parcels after optimisation: ", parcels_benchmark[a,:BS], 
                  " / capacity_used: ", cap_used[a,:BS],  " / time: ",round(time_benchmark[a,:BS], digits = 3))
        end

        ## Start our reproduction of the  bestselling heuristic with an additional local search by
        ## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
        if  start[1,:BSLOC] == 1
            time_benchmark[a,:BSLOC] += @elapsed W = BESTSELLING(trans, Q::Matrix{Int64},
                                                                 capacity::Array{Int64,1})
            time_benchmark[a,:BSLOC] += @elapsed W = LOCALSEARCH(W::Array{Int64,2}, Q::Array{Int64,2})
            cap_used[a,:BSLOC] = sum(W)
            parcels_benchmark[a,:BSLOC] = PARCELSSEND(trans, W, capacity, combination)
            print("\n       bs+loc: parcels after optimisation: ", parcels_benchmark[a,:BSLOC], 
                  " / capacity_used: ", cap_used[a,:BSLOC],  " / time: ",round(time_benchmark[a,:BSLOC], digits = 3))
        end

        ## Start the search for optimal solution with the solver CPLEX
        ## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
        ## FULLOPTUEQ if SKUs can be allocated multiple times
        if start[1,:OPT] == 1
            if sum(capacity) == size(trans,2)
                time_benchmark[a,:OPT] += @elapsed W,gap_optimisation[a,5],popt = FULLOPTEQ(trans::SparseMatrixCSC{Float64, Int64},
                                                                                           capacity::Array{Int64,1},
                                                                                           abort::Int64,
                                                                                           show_opt::Int64,
                                                                                           cpu_cores::Int64,
                                                                                           allowed_gap::Float64,
                                                                                           max_nodes::Int64)
            else
                time_benchmark[a,:OPT] += @elapsed W,gap_optimisation[a,5],popt = FULLOPTUEQ(trans::SparseMatrixCSC{Float64, Int64},
                                                                                            capacity::Array{Int64,1},
                                                                                            abort::Int64,
                                                                                            show_opt::Int64,
                                                                                            cpu_cores::Int64,
                                                                                            allowed_gap::Float64,
                                                                                            max_nodes::Int64)
            end
            cap_used[a,:OPT] = sum(W)
            parcels_benchmark[a,:OPT] = PARCELSSEND(trans, W, capacity, combination)
            print("\n Solution of optimisation: ",popt, "/ Solution of parcel simulation: ", sum(parcel),
            " / time: ", round(time_benchmark[a,:OPT], digits = 3))
        end

        ## Benchmark the random allocation of SKUs
        time_benchmark[a,:RND] += @elapsed parcels_benchmark[a,:RND] = RANDOMBENCH(trans::SparseMatrixCSC{Float64, Int64},
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