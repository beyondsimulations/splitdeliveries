using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Read the CSV file
df = CSV.read("results/overall_results.csv", DataFrame)

# Define mapping from mode names to table headers
mode_mapping = Dict(
    "OPT" => "OPT",
    "QMK" => "QMK-OPT",
    "QMKJ" => "QMK",
    "CHI_0.01" => "CHI-NL",
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

# Define heuristics in order for table
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "OPT"]

# Get unique SKU levels and sort them
sku_levels = sort(unique(filtered_df.skus))

# Calculate results for each SKU level and heuristic
results = DataFrame(skus=Int[], heuristic=String[], avg_duration=Any[], success_rate=Any[])

println("Computation time analysis:")
println("="^50)

for sku in sku_levels
    println("\nSKU level: $sku")

    for heur in heuristics
        subset = filtered_df[(filtered_df.skus.==sku).&(filtered_df.heuristic.==heur), :]

        if nrow(subset) > 0
            # Calculate success rate and average duration
            successful_runs = subset.duration[subset.duration.<3600]  # Assuming 3600s timeout
            total_runs = nrow(subset)
            success_rate = length(successful_runs) / total_runs * 100

            if length(successful_runs) > 0
                avg_duration = round(mean(successful_runs), digits=1)
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

        push!(results, (skus=sku, heuristic=heur, avg_duration=avg_duration, success_rate=success_rate_str))
    end
end

# Display the table
println("\n\n")
println("\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrr}")
println("\\toprule")
print("\$|\\mathcal{I}|\$")
for heur in heuristics
    print(" & $heur")
end
println(" \\\\")
println("\\midrule")

# Print data rows
for sku in sku_levels
    sku_results = results[results.skus.==sku, :]
    if nrow(sku_results) > 0
        # Computation times row
        print("$sku")
        for heur in heuristics
            heur_data = sku_results[sku_results.heuristic.==heur, :]
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
        print(" ")
        for heur in heuristics
            heur_data = sku_results[sku_results.heuristic.==heur, :]
            if nrow(heur_data) > 0
                success_rate = heur_data[1, :success_rate]
                if success_rate == ""
                    print(" & ")
                else
                    print(" & {\\footnotesize \$$success_rate\$}")
                end
            else
                print(" & ")
            end
        end
        println(" \\\\")
        println("\\midrule")
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
