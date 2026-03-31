function FULLOPTUEQ(
    trans::SparseMatrixCSC{Bool,Int64},
    capacity::Array{Int64,1},
    abort::Int64,
    show_opt::Bool,
    cpu_cores::Int64,
    allowed_gap::Float64,
    max_nodes::Int64,
)
    N = size(trans, 2)
    M = size(trans, 1)
    W = count(x -> x > 0, capacity)
    if Int128(M) * N * W > 500_000_000
        error(
            "FULLOPTUEQ: Problem size M×N×W = $(M)×$(N)×$(W) = $(Int128(M)*N*W) exceeds feasible limit for Z variables.",
        )
    end

    # Convert sparse matrix to matrix
    trans = Matrix(trans)

    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity; rev = true)
    allocation = Model(Gurobi.Optimizer)

    if show_opt == false
        set_silent(allocation)
    end

    set_optimizer_attribute(allocation, "MIPGap", allowed_gap)
    set_optimizer_attribute(allocation, "TimeLimit", abort)
    set_optimizer_attribute(allocation, "Threads", cpu_cores)
    set_optimizer_attribute(allocation, "NodeLimit", max_nodes)

    GJ = 1:size(trans, 1)
    GI = 1:size(trans, 2)
    GK = 1:size(capacity, 1)
    @variable(allocation, X[GI, GK], Bin)    ## sku-warehouse allocation
    @variable(allocation, Y[GJ, GK], Bin)    ## send-parcels
    @variable(allocation, Z[GJ, GI, GK], Bin) ## transactions-parcel-allocation
    @objective(allocation, Min, sum(Y[j, k] for j in GJ, k in GK))
    @constraint(
        allocation,
        trans_to_sku[j in GJ, i in GI],
        sum(Z[j, i, k] for k in GK) == trans[j, i]
    )
    @constraint(
        allocation, parcel_to_trans[i in GI, j in GJ, k in GK], Z[j, i, k] <= Y[j, k]
    )
    @constraint(
        allocation, sku_to_warehouse[i in GI, j in GJ, k in GK], Z[j, i, k] <= X[i, k]
    )
    @constraint(allocation, capconstraint[k in GK], sum(X[i, k] for i in GI) == capacity[k])
    @constraint(allocation, minallocation[i in GI], sum(X[i, k] for k in GK) >= 1)
    JuMP.optimize!(allocation)
    G =
        abs(
            objective_bound(allocation)-objective_value(allocation)
        )/abs(objective_value(allocation)+0.00000000001)
    P = objective_value(allocation)
    out = Array{Bool,2}(undef, size(trans, 2), size(capacity, 1)) .= 0
    for i in axes(out, 1)
        for j in axes(out, 2)
            if value.(X[i, j]) > 0
                out[i, j] = 1
            end
        end
    end

    # Check if the solution is empty (no allocations made)
    if sum(out) == 0
        if termination_status(mqkp) == MOI.INFEASIBLE
            error(
                "MQKP optimization resulted in an infeasible model. Please check input parameters.",
            )
        else
            println(
                "MQKP optimization resulted in an empty solution. Status: $(termination_status(mqkp))",
            )
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
                error(
                    "Cannot create a feasible allocation with the given capacities and weights.",
                )
            end
        end
    end

    return out::Array{Bool,2}, G::Float64, P::Float64
end
