## Shared utility functions for Lin et al. (2025) heuristics
## "Multi-Warehouse Assortment Selection: Minimizing Order Splitting in E-Commerce Logistics"

## MCIRANKING: Rank products by decreasing marginal choice probability (order frequency)
## Paper formula: ω_n = Σ π(T) · I(n ∈ T) = sum(trans[:,n]) / num_orders
## Returns a vector of product indices sorted by decreasing omega_n
function MCIRANKING(trans::SparseMatrixCSC{Bool,Int64})
    N = size(trans, 2)
    num_orders = size(trans, 1)

    # Compute marginal choice probability for each product
    omega = zeros(Float64, N)
    for n in 1:N
        omega[n] = sum(trans[:, n]) / num_orders
    end

    # Return product indices sorted by decreasing omega
    return sortperm(omega; rev = true)
end

## OFRM_SUBPROBLEM: Solve the Order Fulfillment Rate Maximization subproblem (Problem 6 in Lin et al.)
## Selects K products for a warehouse to maximize the number of orders fulfillable,
## given that some products must be included (those not in the other warehouse).
## IMPORTANT: K is a cardinality count (number of products). This formulation assumes
## uniform unit weights. Callers must ensure sku_weight is uniform before invoking.
##
## Variables: x[1:N] binary (product selection), y[1:M_uncovered] continuous (order fulfillment)
## Objective: max sum(y_j) -- maximize fulfilled uncovered orders
## Constraints: sum(x) = K, y_j <= x_i for each product i in order j, x_i=1 for must_include
function OFRM_SUBPROBLEM(
    trans::SparseMatrixCSC{Bool,Int64},
    uncovered_orders::Vector{Int64},
    K::Int64,
    must_include::Vector{Int64},
    candidate_products::Vector{Int64},
    solver_name::String,
    abort::Int64,
    show_opt::Bool,
    cpu_cores::Int64,
    allowed_gap::Float64,
    max_nodes::Int64,
)
    N_candidates = length(candidate_products)
    M_uncovered = length(uncovered_orders)

    # Defensive check: must_include cannot exceed capacity
    if length(must_include) > K
        error(
            "OFRM_SUBPROBLEM: |must_include| = $(length(must_include)) > K = $K. Infeasible.",
        )
    end

    # Short-circuit: if no uncovered orders, just pick K most popular products
    if M_uncovered == 0
        selected = copy(must_include)
        # Fill remaining capacity with candidate products not already included
        remaining = setdiff(candidate_products, must_include)
        for p in remaining
            if length(selected) >= K
                break
            end
            push!(selected, p)
        end
        return selected, 0.0
    end

    # Build solver model
    if solver_name == "Gurobi"
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "MIPGap", allowed_gap)
        set_optimizer_attribute(model, "TimeLimit", abort)
        set_optimizer_attribute(model, "Threads", cpu_cores)
        set_optimizer_attribute(model, "NodeLimit", max_nodes)
    elseif solver_name == "scip"
        model = Model(SCIP.Optimizer)
        set_optimizer_attribute(model, "limits/gap", allowed_gap)
        set_optimizer_attribute(model, "limits/time", abort)
        set_optimizer_attribute(model, "parallel/maxnthreads", cpu_cores)
    else
        error(
            "OFRM_SUBPROBLEM: Only Gurobi or SCIP solvers are supported. Got: $solver_name"
        )
    end

    if !show_opt
        set_silent(model)
    end

    # Create index mapping: candidate_products[idx] -> product index
    # x[idx] = 1 means candidate_products[idx] is selected
    prod_to_idx = Dict{Int64,Int64}()
    for (idx, p) in enumerate(candidate_products)
        prod_to_idx[p] = idx
    end

    # Variables
    @variable(model, x[1:N_candidates], Bin)
    @variable(model, 0 <= y[1:M_uncovered] <= 1)

    # Capacity constraint: exactly K products selected
    @constraint(model, sum(x) == K)

    # Must-include constraints
    for p in must_include
        if haskey(prod_to_idx, p)
            @constraint(model, x[prod_to_idx[p]] == 1)
        end
    end

    # Order fulfillment constraints: y_j <= x_i for each product i in uncovered order j
    # An order j can only be fulfilled (y_j=1) if ALL its products are selected
    for (j_idx, order_idx) in enumerate(uncovered_orders)
        # Get all products in this order that are candidates
        row = trans[order_idx, :]
        products_in_order = findnz(row)[1]
        for p in products_in_order
            if haskey(prod_to_idx, p)
                @constraint(model, y[j_idx] <= x[prod_to_idx[p]])
            else
                # Product not a candidate (not in this warehouse's scope) -> order cannot be fulfilled
                @constraint(model, y[j_idx] == 0)
                break
            end
        end
    end

    # Objective: maximize fulfilled orders
    @objective(model, Max, sum(y))

    JuMP.optimize!(model)

    # Extract gap
    gap = 0.0
    try
        gap =
            abs(objective_bound(model) - objective_value(model)) /
            abs(objective_value(model) + 1e-11)
    catch
        gap = 0.0
    end

    # Extract selected products
    selected = Int64[]
    for (idx, p) in enumerate(candidate_products)
        if value(x[idx]) > 0.5
            push!(selected, p)
        end
    end

    # Free solver model to prevent memory accumulation in IIH loop
    empty!(model)

    return selected, gap
end
