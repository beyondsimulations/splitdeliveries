# K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) doi:10.1016/j.ejor.2019.07.004
function KLINKS(trans::SparseMatrixCSC{Float64, Int64},
                capacity::Array{Int64,1},
                trials::Int64,
                stagnant::Int64,
                strategy::Int64,
                klinkstatus::Int64)
    trans = Matrix(trans)         
    if CHECKCAPACITY(trans,capacity) == 1
        ## Calculation of links based on orders
        ov = LINKADJUST(trans::Array{Int64,2})
        L = LINKS(trans::Array{Int64,2},ov::Array{Float64,1})
        ## Set trial = 0
        besttrial = 0
        lw_trial = Array{Float64,1}(undef,trials) .= 0
        cw_trial = Array{Int64,3}(undef,size(trans,2),size(capacity,1),trials) .= 0
        if klinkstatus == 1
            print("\nstarting k-link trials: ")
        end
        for y = 1:trials
            ## Random initialisation of the category distribution
            X = RANDOMALLOCONCE(trans::Array{Int64,2},capacity::Array{Int64,1})
            ## Get the matrix CW trial, Calculate the LW trial
            lw_trial[y] = LW(L::Array{Float64,2},X::Array{Int64,2})
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
                    X = STRATEGY1(X::Array{Int64,2},
                                  m::Int64,
                                  L::Array{Float64,2},
                                  capacity::Array{Int64,1},
                                  stop::Int64)
                else
                    X = STRATEGY2(X::Array{Int64,2},
                                  m::Int64,
                                  L::Array{Float64,2},
                                  capacity::Array{Int64,1},
                                  stop::Int64)
                end
            end
            # Save the results from the trial
            lw_trial[y] = LW(L::Array{Float64,2},X::Array{Int64,2})
            cw_trial[:,:,y] = X
        end
        # Save the best result from all trails for export
        besttrial     = findmax(lw_trial)[2]
        cw_best       = cw_trial[:,:,besttrial]
        return cw_best
    end
end
