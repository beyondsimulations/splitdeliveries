function MQKP(trans::SparseMatrixCSC{Bool,Int64},
    capacity::Array{Int64,1},
    sku_weight::Vector{<:Real},
    abort::Int64,
    solv::String,
    show_opt::Bool,
    cpu_cores::Int64,
    allowed_gap::Float64,
    max_nodes::Int64,
    mode::String)

    # Create the Matrix for the objective of the optimisation
    if mode == "QMK"
        Q = COAPPEARENCE(trans, sku_weight)
    elseif mode == "KLINK"
        Q = LINKS(trans, LINKADJUST(trans)) .* (1)
    else
        error("Sorry, only the mode QMK or KLINK is allowed.")
    end

    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)

    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)

    if solv == "Gurobi"
        # Start building the model
        mqkp = Model(Gurobi.Optimizer)

        # Gurobi-specific parameter settings
        set_optimizer_attribute(mqkp, "MIPGap", allowed_gap)
        set_optimizer_attribute(mqkp, "TimeLimit", abort)
        set_optimizer_attribute(mqkp, "Threads", cpu_cores)
        set_optimizer_attribute(mqkp, "NodeLimit", max_nodes)

    elseif solv == "Juniper"
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false, "threads" => cpu_cores)
        mqkp = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                "nl_solver" => ipopt,
                "mip_solver" => highs,
            ),
        )
        set_optimizer_attribute(mqkp, "mip_gap", allowed_gap)
        set_optimizer_attribute(mqkp, "time_limit", abort)

    elseif solv == "scip"
        # Start building the model
        mqkp = Model(SCIP.Optimizer)

        set_optimizer_attribute(mqkp, "limits/gap", allowed_gap)
        set_optimizer_attribute(mqkp, "limits/time", abort)
        set_optimizer_attribute(mqkp, "parallel/maxnthreads", cpu_cores)
    else
        error("Sorry, only the solver SCIP, Gurobi or Juniper is allowed.")
    end

    if show_opt == false
        set_silent(mqkp)
    end

    GI = 1:size(Q, 2)
    GK = 1:length(capacity)
    @variable(mqkp, X[GI, GK], Bin)
    @objective(mqkp, Max, sum(X[i, k] .* X[j, k] .* Q[i, j] for i in GI, j in 1:i-1, k in GK))
    @constraint(mqkp, capconstraint[k in GK], sum(X[i, k] * sku_weight[i] for i in GI) == capacity[k])
    if mode == "QMK"
        @constraint(mqkp, minallocation[i in GI], sum(X[i, k] for k in GK) >= 1)
    elseif mode == "KLINK"
        @constraint(mqkp, minallocation[i in GI], sum(X[i, k] for k in GK) == 1)
    end
    JuMP.optimize!(mqkp)
    G = 0
    try
        G = abs(objective_bound(mqkp) - objective_value(mqkp)) / abs(objective_value(mqkp) + 0.00000000001)
    catch
        G = 0
    end
    out = Array{Bool,2}(undef, size(Q, 2), size(capacity, 1)) .= 0
    for i in axes(out, 1)
        for j in axes(out, 2)
            if value.(X[i, j]) > 0.1
                out[i, j] = 1
            end
        end
    end

    # Check if the solution is empty (no allocations made)
    if sum(out) == 0
        if termination_status(mqkp) == MOI.INFEASIBLE
            error("MQKP optimization resulted in an infeasible model. Please check input parameters.")
        else
            println("MQKP optimization resulted in an empty solution. Status: $(termination_status(mqkp))")
            println("Attempting basic allocation as fallback...")
            # Simple greedy allocation - assign each item to the first warehouse with capacity
            remaining_capacity = copy(capacity)
            for i in GI
                for k in GK
                    if remaining_capacity[k] >= sku_weight[i]
                        out[i, k] = 1
                        remaining_capacity[k] -= sku_weight[i]
                        break
                    end
                end
            end
            # If we still have an empty solution, it's truly infeasible
            if sum(out) == 0
                error("Cannot create a feasible allocation with the given capacities and weights.")
            end
        end
    end

    return out, G
end
