## Iterative Improvement Heuristic (IIH) by Lin et al. (2025)
## "Multi-Warehouse Assortment Selection: Minimizing Order Splitting in E-Commerce Logistics"
## Implements Algorithm 1 from the main paper (2-warehouse overlapping only)

function IIH(
    trans::SparseMatrixCSC{Bool,Int64},
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real},
    abort::Int64,
    solver_name::String,
    show_opt::Bool,
    cpu_cores::Int64,
    allowed_gap::Float64,
    max_nodes::Int64,
    max_iih_iterations::Int64,
    epsilon_iih::Float64,
)

    # IIH requires uniform sku_weight (unit weights only)
    if !all(w -> w == sku_weight[1], sku_weight)
        error("IIH: Requires uniform sku_weight. Got non-uniform weights.")
    end

    # Sort warehouses by decreasing capacity
    capacity = sort(capacity; rev = true)
    D = length(capacity)
    N = size(trans, 2)
    M = size(trans, 1)

    if D != 2
        error("IIH: Only 2-warehouse configurations are supported. Got $D warehouses.")
    end

    K1 = capacity[1]
    K2 = capacity[2]

    # Initialize with eMCI allocation
    W = EMCIALLOC(trans, capacity, sku_weight)

    # --- Initial step (Paper Algorithm 1, Step 1): Fix S1, optimize S2 via OFRM ---
    S1_init = findall(W[:, 1] .== true)
    must_include_2_init = setdiff(1:N, S1_init)
    candidate_2_init = collect(1:N)

    uncovered_by_1_init = Int64[]
    for j in 1:M
        order_products = findnz(trans[j, :])[1]
        if !all(p -> W[p, 1], order_products)
            push!(uncovered_by_1_init, j)
        end
    end

    selected_2_init, _ = OFRM_SUBPROBLEM(
        trans,
        uncovered_by_1_init,
        K2,
        collect(must_include_2_init),
        candidate_2_init,
        solver_name,
        abort,
        show_opt,
        cpu_cores,
        allowed_gap,
        max_nodes,
    )

    W[:, 2] .= false
    for p in selected_2_init
        W[p, 2] = true
    end

    uncovered_by_1_init = nothing
    selected_2_init = nothing

    # Compute initial split deliveries using warehouse combinations
    combination = COMBINEWAREHOUSES(capacity)
    best_splits = PARCELSSEND(trans, W, capacity, combination)
    total_gap = 0.0
    iteration = 0

    for iter in 1:max_iih_iterations
        improved = false

        # --- Step 1: Fix S2, re-optimize S1 ---
        S2 = findall(W[:, 2] .== true)  # products in warehouse 2
        S1_old = findall(W[:, 1] .== true)

        # Products NOT in S2 must be in S1 (every product needs at least one warehouse)
        must_include_1 = setdiff(1:N, S2)
        # Products in S2 are candidates for also being in S1
        candidate_1 = collect(1:N)

        # Orders already fully covered by warehouse 2 alone
        covered_by_2 = Int64[]
        uncovered_by_2 = Int64[]
        for j in 1:M
            order_products = findnz(trans[j, :])[1]
            if all(p -> W[p, 2], order_products)
                push!(covered_by_2, j)
            else
                push!(uncovered_by_2, j)
            end
        end

        selected_1, gap_1 = OFRM_SUBPROBLEM(
            trans,
            uncovered_by_2,
            K1,
            collect(must_include_1),
            candidate_1,
            solver_name,
            abort,
            show_opt,
            cpu_cores,
            allowed_gap,
            max_nodes,
        )

        # Update W for warehouse 1
        W[:, 1] .= false
        for p in selected_1
            W[p, 1] = true
        end

        # Free intermediates from step 1
        covered_by_2 = nothing
        uncovered_by_2 = nothing
        selected_1 = nothing

        # --- Step 2: Fix S1, re-optimize S2 ---
        S1 = findall(W[:, 1] .== true)

        # Products NOT in S1 must be in S2
        must_include_2 = setdiff(1:N, S1)
        candidate_2 = collect(1:N)

        # Orders already fully covered by warehouse 1 alone
        uncovered_by_1 = Int64[]
        for j in 1:M
            order_products = findnz(trans[j, :])[1]
            if !all(p -> W[p, 1], order_products)
                push!(uncovered_by_1, j)
            end
        end

        selected_2, gap_2 = OFRM_SUBPROBLEM(
            trans,
            uncovered_by_1,
            K2,
            collect(must_include_2),
            candidate_2,
            solver_name,
            abort,
            show_opt,
            cpu_cores,
            allowed_gap,
            max_nodes,
        )

        # Update W for warehouse 2
        W[:, 2] .= false
        for p in selected_2
            W[p, 2] = true
        end

        # Free intermediates from step 2
        uncovered_by_1 = nothing
        selected_2 = nothing

        total_gap = max(gap_1, gap_2)

        # Evaluate improvement
        current_splits = PARCELSSEND(trans, W, capacity, combination)
        improvement = best_splits - current_splits
        iteration = iter

        if improvement > epsilon_iih
            best_splits = current_splits
            improved = true
        end

        GC.gc(false)

        if !improved
            break
        end
    end

    combination = nothing

    return W, total_gap, iteration
end
