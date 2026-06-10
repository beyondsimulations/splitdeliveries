## Extended MCI heuristic for D warehouses by Lin et al. (2025)
## "Multi-Warehouse Assortment Selection: Minimizing Order Splitting in E-Commerce Logistics"
## Uniform weights: implements Algorithm 5 from Appendix EC.5 exactly (closed-form replication)
## Non-uniform weights: greedy weight-aware replication ranked by MCI / sku_weight

function EMCIALLOC(
    trans::SparseMatrixCSC{Bool,Int64}, capacity::Vector{<:Real}, sku_weight::Vector{<:Real}
)
    D = length(capacity)
    N = size(trans, 2)
    uniform_weights = all(w -> w == sku_weight[1], sku_weight)

    W = zeros(Bool, N, D)

    if uniform_weights
        # --- Closed-form Algorithm 5 (original EMCI) ---
        capacity = sort(Int64.(capacity); rev = true)
        ranking = MCIRANKING(trans)
        N_tilde = sum(capacity) - N

        if N_tilde <= 0
            capacity_left = copy(capacity)
            for pos in 1:N
                product = ranking[pos]
                for d in 1:D
                    if capacity_left[d] > 0
                        W[product, d] = true
                        capacity_left[d] -= 1
                        break
                    end
                end
            end
            return W
        end

        # Replication levels via closed-form formula
        r = zeros(Int64, D - 1)
        for j in 1:(D - 1)
            numerator = N_tilde - sum((D - i) * r[i] for i in 1:(j - 1); init = 0)
            denominator = D - j
            r[j] = max(0, floor(Int64, numerator / denominator))
        end

        # Phase 1: Assign replicated products in MCI order
        pos = 1
        for j in 1:(D - 1)
            num_warehouses = D - j + 1
            for _ in 1:r[j]
                if pos > N
                    break
                end
                product = ranking[pos]
                for d in 1:num_warehouses
                    W[product, d] = true
                end
                pos += 1
            end
        end

        # Phase 2: Assign remaining products to single warehouses
        capacity_used = vec(sum(W; dims = 1))
        capacity_left = capacity .- capacity_used
        for d in 1:D
            slots = capacity_left[d]
            for _ in 1:slots
                if pos > N
                    break
                end
                product = ranking[pos]
                W[product, d] = true
                pos += 1
            end
        end

    else
        # --- Greedy weight-aware replication ---
        capacity = sort(Float64.(capacity); rev = true)

        # Rank by MCI / weight (benefit per unit capacity consumed)
        num_orders = size(trans, 1)
        mci = zeros(Float64, N)
        for n in 1:N
            mci[n] = sum(trans[:, n]) / num_orders
        end
        ranking = sortperm([mci[i] / sku_weight[i] for i in 1:N]; rev = true)

        budget = sum(capacity) - sum(sku_weight)
        if budget < -1e-6
            error(
                "EMCIALLOC: Insufficient capacity ($(sum(capacity))) for total weight ($(sum(sku_weight))).",
            )
        end
        budget = max(0.0, budget)

        cap_left = copy(capacity)

        # Per-warehouse replication budget
        repl_budget = budget / (D - 1)

        # Phase 1: Replicate top-ranked SKUs to ALL warehouses within budget
        pos = 1
        budget_left = repl_budget
        while pos <= N && budget_left > 0
            product = ranking[pos]
            w = sku_weight[product]
            if w > budget_left
                break
            end
            # Check all warehouses can fit this SKU
            if !all(d -> cap_left[d] >= w, 1:D)
                break
            end
            # Check feasibility before committing
            for d in 1:D
                cap_left[d] -= w
            end
            remaining_w = [
                sku_weight[ranking[idx]] for
                idx in (pos + 1):N if sum(W[ranking[idx], :]) == 0
            ]
            if !FEASIBLE_REMAINING(cap_left, remaining_w)
                for d in 1:D
                    cap_left[d] += w
                end
                break
            end
            # Commit replication
            for d in 1:D
                W[product, d] = true
            end
            budget_left -= w
            pos += 1
        end

        # Phase 2: Assign remaining unallocated SKUs to single best-fit warehouse
        for idx in pos:N
            product = ranking[idx]
            if sum(W[product, :]) > 0
                continue
            end
            w = sku_weight[product]
            d_best = argmax(cap_left)
            if cap_left[d_best] < w
                SWAPREPAIR!(W, cap_left, sku_weight, product)
            else
                W[product, d_best] = true
                cap_left[d_best] -= w
            end
        end
    end

    # Verify all products allocated
    for n in 1:N
        if sum(W[n, :]) < 1
            error("EMCIALLOC: Product $n not allocated to any warehouse.")
        end
    end

    return W
end
