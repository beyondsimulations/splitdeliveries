function MQKP(trans::Array{Int64,2},
              capacity::Array{Int64,1},
              Q::Array{Int64,2},
              abort::Int64,
              solv::String,
              show_opt::Int64,
              cpu_cores::Int64,
              allowed_gap::Float64,
              max_nodes::Int64)

    if CHECKCAPACITY(trans,capacity) == 1
        mqkp = Model(GAMS.Optimizer)
        if show_opt == 0
            set_silent(mqkp)
        end
        set_optimizer_attribute(mqkp, GAMS.ModelType(), "MIQCP")
        set_optimizer_attribute(mqkp, "Solver",  solv)
        set_optimizer_attribute(mqkp, "OptCR",   allowed_gap)
        set_optimizer_attribute(mqkp, "ResLim",  abort)
        set_optimizer_attribute(mqkp, "Threads", cpu_cores)
        set_optimizer_attribute(mqkp, "NodLim",  max_nodes)
        GI = 1:size(trans,2)
        GK = 1:size(capacity,1)
        @variable(mqkp, X[GI,GK], Bin)
        @objective(mqkp, Max, sum(X[i,k] .* X[j,k] .* Q[i,j] for i in GI, j in 1:i-1, k in GK))
        @constraint(mqkp, capconstraint[k in GK], sum(X[i,k] for i in GI) == capacity[k])
        @constraint(mqkp, minallocation[i in GI], sum(X[i,k] for k in GK) >= 1)
        JuMP.optimize!(mqkp)
        G = abs(objective_bound(mqkp)-objective_value(mqkp))/abs(objective_value(mqkp)+0.00000000001)
        out = Array{Int64,2}(undef,size(trans,2),size(capacity,1)) .= 0
        for i = 1:size(out,1)
            for j = 1:size(out,2)
                if value.(X[i,j]) > 0
                    out[i,j] = 1
                end
            end
        end
        return out::Array{Int64,2},G::Float64
    end
end