## import packages
include("load_packages.jl")

experiment = "3_1000skus"
capa = 1000:1000:15000
ware = 2:1:4
diff = 0.0:0.10:0.20
buff = 0.0:0.10:0.20

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
                    vararray .= round(capacity/wstep)
                else
                    max_w  = round((1+dstep)*capacity/wstep)
                    step_w = - round(((2*dstep)/(wstep-1))*capacity/wstep)
                    min_w  = round((1+dstep)*capacity/wstep) + (wstep-1) * step_w
                    vstep = max_w:step_w:min_w
                    for i = 1:wstep
                        vararray[i] = vstep[i]
                    end
                end
                for warehouse = 1:wstep
                    if warehouse != wstep
                        benchmark[current_row,warehouse] = vararray[warehouse]
                    else
                        benchmark[current_row,warehouse] = round(Int64, capacity - sum(benchmark[current_row,:]))
                    end
                end
                benchmark[current_row,:] = sort(benchmark[current_row,:], rev=true)
                sku_list[current_row] = cstep
                global current_row += 1
            end
        end
    end
end
writedlm("capacity/capacity_$experiment.csv", benchmark, ';')
writedlm("capacity/skus_$experiment.csv", sku_list)

