# script to generate new problem instances for benchmarking
## import packages
include("load_packages.jl")

skus = 20000:20000:100000
experiment = "large"
ware = 2:2:10
diff = 0.0:0.2:0.2
buff = 0.0:0.2:0.2

total_rows = length(skus)*length(ware)*length(diff)*length(buff)
total_cols = maximum(ware)

benchmark = Matrix{Int64}(undef,total_rows,total_cols) .= 0
sku_list  = Vector{Int64}(undef,total_rows) .= 0
differences = Vector{Float64}(undef,total_rows) .= 0
buffers = Vector{Float64}(undef,total_rows) .= 0

current_row = 1
for cstep in skus
    for bstep in buff
        for wstep in ware
            for dstep in diff
                skuscity = cstep * (1+bstep)
                vararray = Vector{Float64}(undef,wstep) .= 0
                if dstep == 0.00
                    vararray .= round(skuscity/wstep)
                else
                    # Calculate average skuscity per warehouse
                    avg_skuscity = skuscity/wstep
                    
                    # First warehouse has diff% more skuscity than average
                    # Last warehouse has diff% less skuscity than average
                    first_warehouse = round(avg_skuscity * (1 + dstep))
                    last_warehouse = round(avg_skuscity * (1 - dstep))
                    
                    # If we have only 2 warehouses
                    if wstep == 2
                        vararray[1] = max(1, first_warehouse)
                        vararray[2] = max(1, last_warehouse)
                    else
                        # Create a linear distribution between first and last
                        step_value = (first_warehouse - last_warehouse) / (wstep - 1)
                        
                        for i = 1:wstep
                            vararray[i] = max(1, round(first_warehouse - (i-1) * step_value))
                        end
                    end
                    
                    # Adjust the total to match required skuscity
                    current_total = sum(vararray)
                    if current_total != skuscity
                        # Distribute the difference across warehouses
                        diff_amount = skuscity - current_total
                        idx = 1
                        while diff_amount != 0
                            adjustment = diff_amount > 0 ? 1 : -1
                            vararray[idx] += adjustment
                            diff_amount -= adjustment
                            idx = (idx % wstep) + 1
                        end
                    end
                end
                assigned_skuscity = 0
                for warehouse = 1:wstep-1
                    benchmark[current_row,warehouse] = vararray[warehouse]
                    assigned_skuscity += vararray[warehouse]
                end
                benchmark[current_row,wstep] = max(0, round(Int64, skuscity - assigned_skuscity))
                benchmark[current_row,:] = sort(benchmark[current_row,:], rev=true)
                sku_list[current_row] = cstep
                differences[current_row] = dstep
                buffers[current_row] = bstep
                global current_row += 1
            end
        end
    end
end
writedlm("capacity/capacity_$experiment.csv", benchmark, ';')
writedlm("capacity/skus_$experiment.csv", sku_list)
writedlm("capacity/diff_$experiment.csv", differences)
writedlm("capacity/buff_$experiment.csv", buffers)

