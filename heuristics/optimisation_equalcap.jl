function FULLOPTEQ(trans::SparseMatrixCSC{Bool, Int64},
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
    GJ = 1:size(trans,1)
    GI = 1:size(trans,2)
    GK = 1:size(capacity,1)
    @variable(allocation, X[GI,GK], Bin)
    @variable(allocation, Y[GJ,GK], Bin)
    @objective(allocation, Min, sum(Y[j,k] for j in GJ, k in GK))
    @constraint(allocation, capconstraint[k in GK], sum(X[i,k] for i in GI) == capacity[k])
    @constraint(allocation, minallocation[i in GI], sum(X[i,k] for k in GK) == 1)
    @constraint(allocation, mapping[j in GJ, k in GK, i in GI; trans[j,i] == 1], X[i,k] <= Y[j,k])
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
    return out::Matrix{Bool},G::Float64,P::Float64
end