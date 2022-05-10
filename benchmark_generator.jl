## import packages
include("load_packages.jl")

experiment = "2_sensitivity"
capa = 400:400:4000
ware = 2:1:3
diff = 0.0:0.05:0.25
buff = 0.0:0.2:0.2

total_rows = length(capa)*length(ware)*length(diff)*length(buff)
total_cols = length(ware)+1

#total_rows = capacity_steps * warehouse_versions * (capacity_diffsteps) * (buffer_addition + 1)
#total_cols = warehouse_versions + 1
#base_cap   = min_cap:(max_cap-min_cap)/(capacity_steps-1):max_cap
#base_ware  = 2:warehouse_versions+1

benchmark = Matrix{Int64}(undef,total_rows,total_cols) .= 0
sku_list  = Vector{Int64}(undef,total_rows) .= 0

current_row = 1
for cstep in capa
    for bstep in buff
        for wstep in ware
            for dstep in diff
                capacity = cstep * (1+bstep)
                vararray = Vector{Float64}(undef,wstep) .= 0
                if dstep == 0.00
                    vararray .= 1/wstep
                else
                    vstep = 1+dstep:-(2*dstep)/wstep:1-dstep
                    for i = 1:wstep
                        vararray[i] = vstep[i]/wstep
                    end
                end
                for warehouse = 1:wstep
                    if warehouse == 1
                        benchmark[current_row,warehouse] = ceil(round(vararray[warehouse] * capacity, digits = 2))
                    elseif  warehouse != wstep
                        benchmark[current_row,warehouse] = round(vararray[warehouse] * capacity)
                    else
                        benchmark[current_row,warehouse] = round(round(capacity - sum(benchmark[current_row,:]), digits = 2))
                    end
                end
                sku_list[current_row] = cstep
                global current_row += 1
            end
        end
    end
end
writedlm("capacity/capacity_$experiment.csv", benchmark, ';')
writedlm("capacity/skus_$experiment.csv", sku_list)

