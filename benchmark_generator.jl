## import packages
include("load_packages.jl")

experiment = "3_calculation_time_test"
capacity_steps      = 2
min_cap             = 25
max_cap             = 50
warehouse_versions  = 3
buffer_addition     = 1
buffer_size         = 0.2
capacity_diffsteps  = 3
capacity_diffmax    = 0.3

total_rows = capacity_steps * warehouse_versions * (capacity_diffsteps) * (buffer_addition + 1)
total_cols = warehouse_versions + 1
base_cap   = min_cap:(max_cap-min_cap)/(capacity_steps-1):max_cap
base_ware  = 2:warehouse_versions+1

benchmark = Matrix{Int64}(undef,total_rows,total_cols) .= 0
sku_list  = Vector{Int64}(undef,total_rows) .= 0

current_row = 1
for buffstep = 1:buffer_addition + 1
    for cstep in base_cap
        if buffstep == 1
            capacity = cstep
        else
            capacity = cstep * (1+buffer_size)
        end
        for warestep in base_ware
            for capvarstep = 1:capacity_diffsteps
                vararray = Vector{Float64}(undef,warestep) .= 0
                if capvarstep == 1
                    vararray .= 1/length(vararray)
                else
                    capacity_diff = capacity_diffmax/(capacity_diffsteps-capvarstep+1)
                    step = 1
                    for i = (1+capacity_diff):-(capacity_diff*2)/(length(vararray)-1):(1-capacity_diff)
                        vararray[step] = i/length(vararray)
                        step += 1
                    end
                end
                for warehouse = 1:warestep
                    if warehouse == 1
                        benchmark[current_row,warehouse] = ceil(vararray[warehouse] * capacity)
                    elseif  warehouse != warestep
                        benchmark[current_row,warehouse] = round(vararray[warehouse] * capacity)
                    else
                        benchmark[current_row,warehouse] = capacity - sum(benchmark[current_row,:])
                    end
                end
                sku_list[current_row] = sum(benchmark[current_row,:])
                global current_row += 1
            end
        end
    end
end
writedlm("capacity/capacity_$experiment.csv", benchmark, ';')
writedlm("capacity/skus_$experiment.csv", sku_list)

