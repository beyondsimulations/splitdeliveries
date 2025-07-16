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
    "BS" => "BS",
    "RND" => "RND"
)

# Filter for relevant modes and map to standard names
filtered_df = df[in.(df.mode, Ref(keys(mode_mapping))), :]
filtered_df.heuristic = [mode_mapping[mode] for mode in filtered_df.mode]



# Add order/SKU ratio and split ratio calculation
filtered_df.order_sku_ratio = filtered_df.orders ./ filtered_df.skus
filtered_df.split_ratio = filtered_df.parcel_test ./ filtered_df.orders

# Define heuristics in order for table (including RND for comparison)
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "RND"]

# Get unique SKU levels and sort them
sku_levels = sort(unique(filtered_df.skus))

# Calculate results for each SKU level, ratio, and heuristic
results = DataFrame(skus=Int[], ratio=Float64[], heuristic=String[], split_ratio=Any[])

println("Average split ratio (in %) depending on the number of SKUs:")
println("="^60)

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
                    # Calculate success rate and average split ratio
                    successful_runs = subset.duration[subset.duration.<3600]
                    successful_split_ratios = subset.split_ratio[subset.duration.<3600]
                    total_runs = nrow(subset)
                    success_rate = length(successful_runs) / total_runs * 100



                    # Only include results if 100% success rate
                    if success_rate == 100.0 && length(successful_split_ratios) > 0
                        # Calculate actual split ratio
                        avg_split_ratio = mean(successful_split_ratios)
                        actual_split_ratio = round(avg_split_ratio * 100, digits=2)


                    else
                        # Skip this result if not 100% success
                        actual_split_ratio = nothing

                    end
                else
                    actual_split_ratio = nothing
                end

                # Only add to results if we have a valid result
                if actual_split_ratio !== nothing
                    push!(results, (skus=sku, ratio=ratio, heuristic=heur, split_ratio=actual_split_ratio))
                end
            end
        end
    end
end

# Display the table
println("\n\n")
println("\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrrrr}")
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
            # Split ratio row
            ratio_results = sku_results[sku_results.ratio.==ratio, :]

            print("$sku & $(round(Int, ratio))")
            for heur in heuristics
                heur_data = ratio_results[ratio_results.heuristic.==heur, :]
                if nrow(heur_data) > 0
                    split_ratio = heur_data[1, :split_ratio]
                    if split_ratio isa Number
                        print(" & \$$split_ratio\$")
                    else
                        print(" & $split_ratio")
                    end
                else
                    print(" & -")
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
println("      \\item \\textit{Notes.} Split ratios as percentages (parcel\\_test/orders × 100) are shown for each algorithm (lower is better). Only algorithms with 100\\% success rates are shown. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.")
println("      \\item Results are only displayed for algorithms that successfully solved all instances within 3,600 seconds.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")
