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
    "CHI_0.01" => "CHI",
    "KL" => "KL",
    "GO" => "GO",
    "GP" => "GP",
    "GS" => "GS",
    "BS" => "BS",
    "RND" => "RND",
)

# Filter for relevant modes and map to standard names
filtered_df = df[in.(df.mode, Ref(keys(mode_mapping))), :]
filtered_df.heuristic = [mode_mapping[mode] for mode in filtered_df.mode]

# Add order/SKU ratio and split ratio calculations for both test and train
filtered_df.order_sku_ratio = filtered_df.orders ./ filtered_df.skus
filtered_df.test_split_ratio = filtered_df.parcel_test ./ filtered_df.orders
filtered_df.train_split_ratio = filtered_df.parcel_train ./ filtered_df.orders

# Define heuristics in order for table (including RND for comparison)
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "RND"]

# Get unique SKU levels and sort them
sku_levels = sort(unique(filtered_df.skus))

# Calculate results for each SKU level, ratio, and heuristic
results = DataFrame(;
    skus = Int[],
    ratio = Float64[],
    heuristic = String[],
    test_split_ratio = Any[],
    train_split_ratio = Any[],
)

println(
    "Average split ratio (in %) depending on the number of SKUs with train/test comparison:"
)
println("="^60)

for sku in sku_levels
    println("\nSKU level: $sku")

    # Get unique ratios for this SKU level
    sku_data = filtered_df[filtered_df.skus .== sku, :]
    if nrow(sku_data) > 0
        unique_ratios = sort(unique(sku_data.order_sku_ratio))
        println("Order/SKU ratios: $(join(unique_ratios, ", "))")

        for ratio in unique_ratios
            for heur in heuristics
                subset = filtered_df[
                    (filtered_df.skus .== sku) .& (filtered_df.heuristic .== heur) .& (filtered_df.order_sku_ratio .== ratio),
                    :,
                ]

                if nrow(subset) > 0
                    # Calculate success rate and average split ratios
                    successful_runs = subset.duration[subset.duration .< 3600]
                    successful_test_ratios = subset.test_split_ratio[subset.duration .< 3600]
                    successful_train_ratios = subset.train_split_ratio[subset.duration .< 3600]
                    total_runs = nrow(subset)
                    success_rate = length(successful_runs) / total_runs * 100

                    # Only include results if 100% success rate
                    if success_rate == 100.0 && length(successful_test_ratios) > 0
                        # Calculate actual split ratios
                        avg_test_ratio = mean(successful_test_ratios)
                        avg_train_ratio = mean(successful_train_ratios)
                        actual_test_ratio = round(avg_test_ratio * 100; digits = 2)
                        actual_train_ratio = round(avg_train_ratio * 100; digits = 2)
                    else
                        # Skip this result if not 100% success
                        actual_test_ratio = nothing
                        actual_train_ratio = nothing
                    end
                else
                    actual_test_ratio = nothing
                    actual_train_ratio = nothing
                end

                # Only add to results if we have a valid result
                if actual_test_ratio !== nothing
                    push!(
                        results,
                        (
                            skus = sku,
                            ratio = ratio,
                            heuristic = heur,
                            test_split_ratio = actual_test_ratio,
                            train_split_ratio = actual_train_ratio,
                        ),
                    )
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
    sku_results = results[results.skus .== sku, :]
    if nrow(sku_results) > 0
        unique_ratios_for_sku = sort(unique(sku_results.ratio))

        for ratio in unique_ratios_for_sku
            # Test split ratio row
            ratio_results = sku_results[sku_results.ratio .== ratio, :]

            print("$sku & $(round(Int, ratio))")
            for heur in heuristics
                heur_data = ratio_results[ratio_results.heuristic .== heur, :]
                if nrow(heur_data) > 0
                    test_ratio = heur_data[1, :test_split_ratio]
                    if test_ratio isa Number
                        print(" & \$$test_ratio\$")
                    else
                        print(" & $test_ratio")
                    end
                else
                    print(" & -")
                end
            end
            println(" \\\\")

            # Train split ratio row
            print(" & {\\footnotesize (train)}")
            for heur in heuristics
                heur_data = ratio_results[ratio_results.heuristic .== heur, :]
                if nrow(heur_data) > 0
                    train_ratio = heur_data[1, :train_split_ratio]
                    if train_ratio isa Number
                        print(" & {\\footnotesize \$$train_ratio\$}")
                    else
                        print(" & {\\footnotesize $train_ratio}")
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
println(
    "      \\item \\textit{Notes.} Split ratios as percentages are shown for each algorithm (lower is better). Test performance is shown in the first row, training performance in the second row (smaller font). Only algorithms with 100\\% success rates are shown. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.",
)
println(
    "      \\item Results are only displayed for algorithms that successfully solved all instances within 3,600 seconds.",
)
println("\\end{tablenotes}")
println("\\end{threeparttable}")
