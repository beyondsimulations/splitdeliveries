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
        Q = COAPPEARENCE(trans,sku_weight)
    elseif mode == "KLINK"
        Q = LINKS(trans,LINKADJUST(trans)).*(1)
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
    else
        error("Sorry, only the solver Gurobi is allowed.")
    end
    if show_opt == false
        set_silent(mqkp)
    end

    GI = 1:size(Q,2)
    GK = 1:length(capacity)
    @variable(mqkp, X[GI,GK], Bin)
    @objective(mqkp, Max, sum(X[i,k] .* X[j,k] .* Q[i,j] for i in GI, j in 1:i-1, k in GK))
    @constraint(mqkp, capconstraint[k in GK], sum(X[i,k] * sku_weight[i] for i in GI) == capacity[k])
    if mode == "QMK"
        @constraint(mqkp, minallocation[i in GI], sum(X[i,k] for k in GK) >= 1)
    elseif mode == "KLINK"
        @constraint(mqkp, minallocation[i in GI], sum(X[i,k] for k in GK) == 1)
    end
    JuMP.optimize!(mqkp)
    G = 0
    try
        G = abs(objective_bound(mqkp)-objective_value(mqkp))/abs(objective_value(mqkp)+0.00000000001)
    catch
        G = 0
    end
    out = Array{Bool,2}(undef,size(Q,2),size(capacity,1)) .= 0
    for i in axes(out,1)
        for j in axes(out,2)
            if value.(X[i,j]) > 0.1
                out[i,j] = 1
            end
        end
    end
    return out,G
end