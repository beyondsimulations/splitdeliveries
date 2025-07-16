using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Read the CSV file
df = CSV.read("results/overall_results.csv", DataFrame)

# Define mapping from mode names to table headers
mode_mapping = Dict(
    "OPT" => "OPT",
    "QMK" => "QMK",
    "QMKJ" => "QMKJ",
    "CHI_0.01" => "CHIL",
    "CHIM_0.01" => "CHI",
    "KL" => "KL",
    "GO" => "GO",
    "GP" => "GP",
    "GS" => "GS",
    "BS" => "BS"
)

# Filter for relevant modes and map to standard names
filtered_df = df[in.(df.mode, Ref(keys(mode_mapping))), :]
filtered_df.heuristic = [mode_mapping[mode] for mode in filtered_df.mode]

# Add order/SKU ratio calculation
filtered_df.order_sku_ratio = filtered_df.orders ./ filtered_df.skus

# Define heuristics in order for table
heuristics = ["QMK", "QMKJ", "CHI", "KL", "GO", "GP", "GS", "BS", "OPT"]

# Get unique SKU levels and sort them
sku_levels = sort(unique(filtered_df.skus))

# Calculate results for each SKU level, ratio, and heuristic
results = DataFrame(skus=Int[], ratio=Float64[], heuristic=String[], avg_duration=Any[], success_rate=Any[])

println("Order/SKU ratio analysis:")
println("="^50)

for sku in sku_levels
    println("\nSKU level: $sku")

    # Get unique ratios for this SKU level
    sku_data = filtered_df[filtered_df.skus.==sku, :]
    if nrow(sku_data) > 0
        unique_ratios = sort(unique(sku_data.order_sku_ratio))
        println("Order/SKU ratios: $(join(unique_ratios, ", "))")

        for ratio in unique_ratios
            for heur in heuristics
                subset = filtered_df[(filtered_df.skus.==sku).&(filtered_df.heuristic.==heur).&(filtered_df.order_sku_ratio.==ratio), :]

                if nrow(subset) > 0
                    # Calculate success rate and average duration
                    successful_runs = subset.duration[subset.duration.<3600]  # Assuming 3600s timeout
                    total_runs = nrow(subset)
                    success_rate = length(successful_runs) / total_runs * 100

                    if length(successful_runs) > 0
                        avg_duration = round(mean(successful_runs), digits=2)
                    else
                        avg_duration = "^a"
                        success_rate = ""  # Empty for ^a cases
                    end

                    # Format success rate
                    if success_rate isa Number
                        if success_rate == 100.0
                            success_rate_str = "100\\%"
                        else
                            success_rate_str = "$(round(Int, success_rate))\\%"
                        end
                    else
                        success_rate_str = success_rate
                    end
                else
                    avg_duration = "^a"
                    success_rate_str = ""
                end

                push!(results, (skus=sku, ratio=ratio, heuristic=heur, avg_duration=avg_duration, success_rate=success_rate_str))
            end
        end
    end
end

# Display the table
println("\n\n")
println("\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrrr}")
println("\\toprule")
print("\$|\\mathcal{I}|\$ & Ratio")
for heur in heuristics
    print(" & $heur")
end
println(" \\\\")
println("\\midrule")

# Print data rows
for sku in sku_levels
    sku_results = results[results.skus.==sku, :]
    if nrow(sku_results) > 0
        unique_ratios_for_sku = sort(unique(sku_results.ratio))

        for ratio in unique_ratios_for_sku
            # Computation times with superscript success rates
            ratio_results = sku_results[sku_results.ratio.==ratio, :]

            print("$sku & $(round(Int, ratio))")
            for heur in heuristics
                heur_data = ratio_results[ratio_results.heuristic.==heur, :]
                if nrow(heur_data) > 0
                    duration = heur_data[1, :avg_duration]

                    if duration == "^a"
                        print(" & \$^a\$")
                    elseif duration isa Number
                        # Format the duration
                        if duration >= 1000
                            formatted = replace(string(duration), r"(?<=\d)(?=(\d{3})+(?!\d))" => ",")
                            print(" & \$$formatted\$")
                        else
                            print(" & \$$duration\$")
                        end
                    else
                        print(" & $duration")
                    end
                else
                    print(" & \$^a\$")
                end
            end
            println(" \\\\")

            # Second row: success rates in smaller font
            print(" & ")
            for heur in heuristics
                heur_data = ratio_results[ratio_results.heuristic.==heur, :]
                if nrow(heur_data) > 0
                    success_rate = heur_data[1, :success_rate]
                    if success_rate == ""
                        print(" & ")
                    else
                        print(" & {\\small \$$success_rate\$}")
                    end
                else
                    print(" & ")
                end
            end
            println(" \\\\")
            println("\\midrule")
        end
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println("      \\item \\textit{Notes.} The computation time is displayed in seconds in the first row, with the success rate shown in the second row. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.")
println("      \\item \$^a\$ Excluded from further benchmarking, as feasible solution could not be computed within 3,600 seconds.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")
