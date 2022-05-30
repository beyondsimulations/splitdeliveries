# K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) doi:10.1016/j.ejor.2019.07.004
function KLINKS(trans::SparseMatrixCSC{Bool, Int64},
                capacity::Array{Int64,1},
                trials::Int64,
                stagnant::Int64,
                strategy::Int64,
                klinkstatus::Int64)      
    if CHECKCAPACITY(trans,capacity) == 1
        ## Calculation of links based on orders
        ov = LINKADJUST(trans::SparseMatrixCSC{Bool, Int64})
        L  = LINKS(trans::SparseMatrixCSC{Bool, Int64}, ov::Array{Float64,1})
        ## Set trial = 0
        besttrial = 0
        lw_trial = Array{Float64,1}(undef,trials) .= 0
        cw_trial = Array{Bool,3}(undef,size(trans,2),size(capacity,1),trials) .= 0
        if klinkstatus == 1
            print("\nstarting k-link trials: ")
        end
        time_limit = Minute(abort/60)
        start = Dates.now()
        time_elapsed = Minute(0)
        for y = 1:trials
            ## Random initialisation of the category distribution
            X = RANDOMALLOCONCE(trans::SparseMatrixCSC{Bool, Int64},capacity::Array{Int64,1})
            ## Get the matrix CW trial, Calculate the LW trial
            lw_trial[y] = LW(L::Array{Float64,2},X::Array{Bool,2})
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
                if strategy == 1
                    X = STRATEGY1(X::Array{Bool,2},
                                  m::Int64,
                                  L::Array{Float64,2},
                                  capacity::Array{Int64,1},
                                  stop::Int64)
                else
                    X = STRATEGY2(X::Array{Bool,2},
                                  m::Int64,
                                  L::Array{Float64,2},
                                  capacity::Array{Int64,1},
                                  stop::Int64)
                end
            end
            # Save the results from the trial
            lw_trial[y] = LW(L::Array{Float64,2},X::Array{Bool,2})
            cw_trial[:,:,y] = X
            time_elapsed = Dates.now() - start
            if time_elapsed > Minute(time_limit)
                print("\nHeuristic stopped due to time elapsed > $time_limit")
                break
            end
        end
        # Save the best result from all trails for export
        besttrial     = findmax(lw_trial)[2]
        cw_best       = cw_trial[:,:,besttrial]
        cw_best = convert(Matrix{Int64},cw_best)
        return cw_best
    end
end
