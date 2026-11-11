# K-LINK heuristic by S. Zhu, X. Hu, K. Huang and Y. Yuan (2021) doi:10.1016/j.ejor.2020.08.024
function KLINKS(
    trans::SparseMatrixCSC{Bool,Int64},
    capacity::Array{Int64,1},
    trials::Int64,
    stagnant::Int64,
    strategy::Int64,
    abort::Int64,
    klinkstatus::Bool,
)
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity; rev = true)
    # Introduce the equal weights for all SKUs
    sku_weight = zeros(Int64, size(trans, 2)) .= 1
    ## Calculation of links based on orders
    N = size(trans, 2)
    if Int128(N) * N > 2_000_000_000
        error(
            "KLINKS: N=$(N) too large for dense N×N LINKS matrix (would require $(round(Int128(N)*N*8/1e9, digits=1)) GB).",
        )
    end
    L = LINKS(trans, LINKADJUST(trans))
    ## Set trial = 0
    besttrial = 0
    lw_trial = zeros(Float64, trials)
    cw_best = zeros(Bool, N, size(capacity, 1))
    if klinkstatus == true
        print("\nstarting k-link trials: ")
    end
    time_limit = Minute(abort/60)
    start = Dates.now()
    time_elapsed = Minute(0)
    local_search_trials = 0
    strategy in (1, 2, 3) || error("Sorry, this strategy is not defined.")
    for y in 1:trials
        ## Random initialisation of the category distribution
        X = RANDOMALLOCONCE(capacity, sku_weight)
        if klinkstatus == true
            print(y, "/ ")
        end
        ## Start the heuristic (Fig. 6 of Zhu et al. 2021): draw a warehouse k, apply the best
        ## improving move of a SKU in k, stop when no warehouse has an improving move left
        # S[i, k]: links of SKU i to the SKUs in warehouse k, so DL(i, m, g) = S[i, g] - S[i, m]
        S = L * X
        wh = [findfirst(view(X, i, :)) for i in 1:N]
        load = vec(sum(X; dims = 1))
        # closed[k]: warehouse k has no improving move in the current allocation; drawing k
        # uniformly among the open warehouses equals the uniform draw of Fig. 6 with the
        # draws of closed warehouses skipped
        closed = falses(size(X, 2))
        # -> 1: relocate the SKU of k with the largest DL into a warehouse with free capacity,
        #       k keeps at least one SKU (Nmin = 1)
        # -> 2: exchange of a SKU of k with a SKU of another warehouse with the largest payoff
        # -> 3: the best relocation of k if it improves, otherwise the best exchange of k
        # a move must gain more than 1e-9 so rounding noise in S cannot make the search cycle
        while !all(closed) && Dates.now() - start <= time_limit
            k = rand(findall(!, closed))
            best = 1e-9
            bi = bj = bg = 0
            if strategy != 2 && load[k] > 1
                for g in axes(S, 2)
                    g != k && load[g] < capacity[g] || continue
                    @inbounds for i in 1:N
                        wh[i] == k || continue
                        dl = S[i, g] - S[i, k]
                        if dl > best
                            best, bi, bg = dl, i, g
                        end
                    end
                end
            end
            if bi == 0 && strategy != 1
                @inbounds for i in 1:N
                    wh[i] == k || continue
                    for j in 1:N
                        g = wh[j]
                        g == k && continue
                        payoff = S[i, g] - S[i, k] + S[j, k] - S[j, g] - 2 * L[j, i]
                        if payoff > best
                            best, bi, bj = payoff, i, j
                        end
                    end
                end
            end
            if bi == 0
                closed[k] = true
            elseif bj == 0
                RELOCATE!(X, S, wh, load, L, bi, bg)
                fill!(closed, false)
            else
                RELOCATE!(X, S, wh, load, L, bi, wh[bj])
                RELOCATE!(X, S, wh, load, L, bj, k)
                fill!(closed, false)
            end
        end
        # Save the results from the trial: LW is half the sum of the links of each SKU
        # within its own warehouse
        lw_trial[y] = sum(S[i, wh[i]] for i in 1:N) / 2
        if besttrial == 0 || lw_trial[y] > lw_trial[besttrial] + 1e-9
            besttrial = y
            copyto!(cw_best, X)
        end
        local_search_trials = y
        time_elapsed = Dates.now() - start
        if time_elapsed > Minute(time_limit)
            print("\nHeuristic stopped due to time elapsed > $time_limit")
            break
        end
        # Stop after stagnant consecutive trials without a better solution
        y - besttrial > stagnant && break
    end
    # Check whether all SKUs are allocated
    if any(y->y < 1, sum(cw_best; dims = 2))
        error("KL: Not all SKUs are allocated!")
    end
    return cw_best, local_search_trials
end
