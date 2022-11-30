# K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) doi:10.1016/j.ejor.2019.07.004
function KLINKS(trans::SparseMatrixCSC{Bool, Int64},
                capacity::Array{Int64,1},
                trials::Int64,
                stagnant::Int64,
                strategy::Int64,
                klinkstatus::Int64)
    # Sort the warehouses by decreasing capacity
    #capacity = sort(capacity, rev=true)   
    if CHECKCAPACITY(trans,capacity) == 1
        ## Calculation of links based on orders
        L  = LINKS(trans,LINKADJUST(trans))
        ## Set trial = 0
        besttrial = 0
        lw_trial = zeros(Float64, trials)
        cw_trial = zeros(Bool,size(trans,2),size(capacity,1),trials)
        if klinkstatus == 1
            print("\nstarting k-link trials: ")
        end
        time_limit = Minute(abort/60)
        start = Dates.now()
        time_elapsed = Minute(0)
        local_search_trials = 0
        for y = 1:trials
            ## Random initialisation of the category distribution
            X = RANDOMALLOCONCE(trans,capacity)
            ## Get the matrix CW trial, Calculate the LW trial
            lw_trial[y] = LW(L,X)
            if klinkstatus == 1
                print(y,"/ ")
            end
            ## Start the heuristic
            stop = 0
            while stop <= stagnant
                # select a random warehouse to test potential improvements
                stop += 1
                m = rand(1:size(X,2))
                # apply either local search strategy 1 or strategy 2
                # -> 1: move SKU to a different warehouse if it improves the objective 
                #       and the warehouse has capacity cap_left
                # -> 2: pair-wise exchange between SKUs if it improves the objective
                # -> 3: both strategies are applied with a chance of 50:50
                if strategy == 1
                    STRATEGY1!(X,m,L,capacity,stop)
                elseif strategy == 2
                    STRATEGY2!(X,m,L,capacity,stop)
                elseif strategy == 3
                    if rand() > 0.5
                        STRATEGY1!(X,m,L,capacity,stop)
                    else
                        STRATEGY2!(X,m,L,capacity,stop)
                    end
                else
                    error("Sorry, this strategy is not defined.")
                end
            end
            # Save the results from the trial
            lw_trial[y] = LW(L,X)
            cw_trial[:,:,y] = X
            local_search_trials = y
            time_elapsed = Dates.now() - start
            if time_elapsed > Minute(time_limit)
                print("\nHeuristic stopped due to time elapsed > $time_limit")
                break
            end
        end
        # Save the best result from all trails for export
        besttrial = findmax(lw_trial)[2]
        cw_best   = cw_trial[:,:,besttrial]
        return cw_best, local_search_trials
    end
end
