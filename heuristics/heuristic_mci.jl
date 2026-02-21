## Extended MCI heuristic for D warehouses by Lin et al. (2025)
## "Multi-Warehouse Assortment Selection: Minimizing Order Splitting in E-Commerce Logistics"
## Implements Algorithm 5 from Appendix EC.5 exactly
## Requires uniform (unit) sku_weight — uses cardinality-based capacity

function EMCIALLOC(trans::SparseMatrixCSC{Bool,Int64},
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real})

    # Algorithm 5 requires uniform sku_weight (cardinality-based capacity)
    if !all(w -> w == sku_weight[1], sku_weight)
        error("EMCIALLOC: Algorithm 5 requires uniform sku_weight. Got non-uniform weights.")
    end

    # Sort warehouses by decreasing capacity
    capacity = sort(Int64.(capacity), rev=true)
    D = length(capacity)
    N = size(trans, 2)

    # Rank products by MCI (decreasing marginal choice probability)
    ranking = MCIRANKING(trans)

    # Compute overlap budget: Ñ = Σ K_d - N (Appendix EC.5)
    N_tilde = sum(capacity) - N

    # Initialize allocation matrix (SKUs x warehouses)
    W = zeros(Bool, N, D)

    if N_tilde <= 0
        # No overlap: partition products across warehouses in MCI order
        # Assign in contiguous blocks by warehouse capacity (largest first)
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

    # Overlapping case: compute replication levels via closed-form formula
    # r[j] = number of products replicated in exactly (D - j + 1) warehouses
    # Formula: r_j = ⌊(Ñ - Σ_{i<j} (D-i)·r_i) / (D-j)⌋
    r = zeros(Int64, D - 1)
    for j in 1:(D-1)
        numerator = N_tilde - sum((D - i) * r[i] for i in 1:(j-1); init=0)
        denominator = D - j
        r[j] = max(0, floor(Int64, numerator / denominator))
    end

    # Phase 1: Assign replicated products in MCI order
    # Products in r[j] go to the (D - j + 1) largest warehouses (warehouses 1..D-j+1)
    pos = 1
    for j in 1:(D-1)
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

    # Phase 2: Assign remaining (non-replicated) products to single warehouses
    # in contiguous blocks by remaining capacity (largest warehouse filled first)
    capacity_used = vec(sum(W, dims=1))
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

    # Verify all products are allocated to at least one warehouse
    for n in 1:N
        if sum(W[n, :]) < 1
            error("EMCIALLOC: Product $n not allocated to any warehouse.")
        end
    end

    return W
end