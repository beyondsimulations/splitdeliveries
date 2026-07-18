using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Read the CSV file
df = CSV.read("results/overall_results.csv", DataFrame)

# Uniform-weight tables only; the weighted modes get dedicated tables.
df = df[df.weight_mode .== "uniform", :]

# Success definition: constructive heuristics count as successful if they
# returned a feasible allocation (every written row is one). For the
# solver-based methods (OPT, QMK) and KL, the run must additionally finish
# within the practical wall-clock cap; the solver time limit is 900 s and
# the cap adds tolerance for model building. OPT further requires a proven
# optimum (gap <= 1e-6).
SUCCESS_CAP = 1.5 * 900
CAPPED_MODES = ["OPT", "QMK", "QMKJ", "KL"]
function is_success(row)
    if row.mode in CAPPED_MODES
        row.duration < SUCCESS_CAP || return false
        row.mode == "OPT" && return row.gap <= 1e-6
        return true
    end
    return true
end

# Define mapping from mode names to table headers
mode_mapping = Dict(
    "OPT" => "OPT",
    "QMK" => "QMK-OPT",
    "QMKJ" => "QMK",
    "CHI_1.0e-5" => "CHI",
    "KL" => "KL",
    "GO" => "GO",
    "GP" => "GP",
    "GS" => "GS",
    "BS" => "BS",
    "EMCI" => "EMCI",
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
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "EMCI", "RND"]

# Get unique SKU levels and sort them
sku_levels = sort(unique(filtered_df.skus))

# Calculate results for each SKU level, ratio, and heuristic
results = DataFrame(;
    skus = Int[],
    ratio = Float64[],
    heuristic = String[],
    test_split_ratio = Any[],
    train_split_ratio = Any[],
    reduced_sample = Bool[],
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

                reduced_sample = false
                if nrow(subset) > 0
                    # Calculate success rate and average split ratios
                    success_mask = is_success.(eachrow(subset))
                    successful_test_ratios = subset.test_split_ratio[success_mask]
                    successful_train_ratios = subset.train_split_ratio[success_mask]
                    total_runs = nrow(subset)
                    success_rate = sum(success_mask) / total_runs * 100

                    if length(successful_test_ratios) > 0
                        # Calculate actual split ratios; mark cells with failures
                        reduced_sample = success_rate < 100.0
                        avg_test_ratio = mean(successful_test_ratios)
                        avg_train_ratio = mean(successful_train_ratios)
                        actual_test_ratio = round(avg_test_ratio * 100; digits = 2)
                        actual_train_ratio = round(avg_train_ratio * 100; digits = 2)
                    else
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
                            reduced_sample = reduced_sample,
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
println("\\begin{tabular}{ll" * "r"^length(heuristics) * "}")
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
                    marker = heur_data[1, :reduced_sample] ? "\$^a\$" : ""
                    if test_ratio isa Number
                        print(" & \$$test_ratio\$$marker")
                    else
                        print(" & $test_ratio$marker")
                    end
                else
                    print(" & -")
                end
            end
            println(" \\\\")

            # Train split ratio row (all heuristics)
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
println(
    "      \\item \\textit{Notes.} Split ratios as percentages are shown for each algorithm (lower is better). Test performance is shown in the first row and training performance in the smaller-font second row. A dash marks methods that are not applicable at the respective scale.",
)
println(
    "      \\item \$^a\$ Mean over successfully solved instances only (reduced sample size, cf.\\ the success definition in Section~5.1).",
)
println("\\end{tablenotes}")
println("\\end{threeparttable}")
