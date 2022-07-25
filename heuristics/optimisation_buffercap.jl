function FULLOPTUEQ(trans::SparseMatrixCSC{Bool, Int64},
                    capacity::Array{Int64,1},
                    abort::Int64,
                    show_opt::Bool,
                    cpu_cores::Int64,
                    allowed_gap::Float64,
                    max_nodes::Int64)
    # Convert sparse matrix to matrix
    trans = Matrix(trans)
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)
    if CHECKCAPACITY(trans,capacity) == 1
        allocation = Model(GAMS.Optimizer)
        if show_opt == 0
            set_silent(allocation)
        end
        set_optimizer_attribute(allocation, GAMS.ModelType(), "MIP")
        set_optimizer_attribute(allocation, "Solver",  "CPLEX")
        set_optimizer_attribute(allocation, "OptCR",   allowed_gap)
        set_optimizer_attribute(allocation, "ResLim",  abort)
        set_optimizer_attribute(allocation, "Threads", cpu_cores)
        set_optimizer_attribute(allocation, "NodLim",  max_nodes)
        GJ in axes(trans,1)
        GI in axes(trans,2)
        GK in axes(capacity,1)
        @variable(allocation, X[GI,GK], Bin)    ## sku-warehouse allocation
        @variable(allocation, Y[GJ,GK], Bin)    ## send-parcels
        @variable(allocation, Z[GJ,GI,GK], Bin) ## transactions-parcel-allocation
        @objective(allocation, Min, sum(Y[j,k] for j in GJ, k in GK))
        @constraint(allocation, trans_to_sku[j in GJ, i in GI], sum(Z[j,i,k] for k in GK) == trans[j,i])
        @constraint(allocation, parcel_to_trans[i in GI, j in GJ, k in GK], Z[j,i,k] <= Y[j,k])
        @constraint(allocation, sku_to_warehouse[i in GI, j in GJ, k in GK], Z[j,i,k] <= X[i,k])
        @constraint(allocation, capconstraint[k in GK], sum(X[i,k] for i in GI) == capacity[k])
        @constraint(allocation, minallocation[i in GI], sum(X[i,k] for k in GK) >= 1)
        JuMP.optimize!(allocation)
        G = abs(objective_bound(allocation)-objective_value(allocation))/abs(objective_value(allocation)+0.00000000001)
        P = objective_value(allocation)
        out = Array{Bool,2}(undef,size(trans,2),size(capacity,1)) .= 0
        for i in axes(out,1)
            for j in axes(out,2)
                if value.(X[i,j]) > 0
                    out[i,j] = 1
                end
            end
        end
        return out::Array{Bool,2},G::Float64,P::Float64
    end
end