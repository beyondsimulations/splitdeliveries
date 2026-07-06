using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Read the CSV file
df = CSV.read("results/overall_results.csv", DataFrame)

# Uniform-weight tables only; the weighted modes get dedicated tables.
df = df[df.weight_mode .== "uniform", :]

# Success definition: a run counts as successful if it returned a feasible
# allocation within the practical wall-clock cap. The solver time limit is
# 900 s; the cap adds tolerance for model building.
SUCCESS_CAP = 1.5 * 900

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
    "RND" => "RND",
)

# Filter for relevant modes and map to standard names
filtered_df = df[in.(df.mode, Ref(keys(mode_mapping))), :]
filtered_df.heuristic = [mode_mapping[mode] for mode in filtered_df.mode]

# Add split ratio calculation
filtered_df.split_ratio = filtered_df.parcel_test ./ filtered_df.orders

# Define key heuristics for comparison (including RND as baseline)
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "RND"]

# Get unique values and sort them
warehouse_levels = sort(unique(filtered_df.wareh))
diff_levels = sort(unique(filtered_df.diff))
buffer_levels = sort(unique(filtered_df.buffer))

println("Complete Warehouse Configuration Analysis")
println("="^60)
println("Warehouses: ", warehouse_levels)
println("Diff levels: ", diff_levels)
println("Buffer levels: ", buffer_levels)

# Calculate results for each combination
results = DataFrame(;
    wareh = Int[],
    diff = Float64[],
    buffer = Float64[],
    heuristic = String[],
    split_ratio = Any[],
    success_rate = Float64[],
    sample_size = Int[],
)

for wareh in warehouse_levels
    for diff in diff_levels
        for buffer in buffer_levels
            for heur in heuristics
                subset = filtered_df[
                    (filtered_df.wareh .== wareh) .& (filtered_df.diff .== diff) .& (filtered_df.buffer .== buffer) .& (filtered_df.heuristic .== heur),
                    :,
                ]

                if nrow(subset) > 0
                    # Calculate success rate and average split ratio
                    successful_runs = subset.duration[subset.duration .< SUCCESS_CAP]
                    successful_ratios = subset.split_ratio[subset.duration .< SUCCESS_CAP]
                    total_runs = nrow(subset)
                    success_rate = length(successful_runs) / total_runs * 100

                    if length(successful_ratios) > 0
                        avg_ratio = mean(successful_ratios)
                        split_ratio = round(avg_ratio * 100; digits = 1)
                        sample_size = length(successful_ratios)
                    else
                        split_ratio = nothing
                        sample_size = 0
                    end

                    push!(
                        results,
                        (
                            wareh = wareh,
                            diff = diff,
                            buffer = buffer,
                            heuristic = heur,
                            split_ratio = split_ratio,
                            success_rate = success_rate,
                            sample_size = sample_size,
                        ),
                    )
                end
            end
        end
    end
end

# Generate the complete table
println("\n\n")
println("\\caption{Split ratio (\\%) by warehouse configuration and algorithm}")
println("\\label{tab:warehouse_complete}")
println("\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrr")
println("\\toprule")
println("Warehouses & Diff & Buffer & QMK & CHI & KL & GO & GP & GS & BS & RND \\\\")
println("\\midrule")

# All four combinations for each warehouse level
config_combinations = [(0.0, 0.0), (0.0, 0.2), (0.2, 0.0), (0.2, 0.2)]

for (w_idx, wareh) in enumerate(warehouse_levels)
    for (c_idx, (diff, buffer)) in enumerate(config_combinations)
        config_results = results[
            (results.wareh .== wareh) .& (results.diff .== diff) .& (results.buffer .== buffer),
            :,
        ]

        # Print warehouse number only for first configuration
        if c_idx == 1
            print("$wareh")
        else
            print(" ")
        end

        # Format diff and buffer values
        diff_str = diff == 0.0 ? "0.0" : string(diff)
        buffer_str = buffer == 0.0 ? "0.0" : string(buffer)
        print(" & $diff_str & $buffer_str")

        for heur in heuristics
            heur_data = config_results[config_results.heuristic .== heur, :]
            if nrow(heur_data) > 0
                ratio = heur_data[1, :split_ratio]
                success_rate = heur_data[1, :success_rate]

                if ratio !== nothing
                    if success_rate < 100.0
                        print(" & $ratio\$^a\$")
                    else
                        print(" & $ratio")
                    end
                else
                    print(" & --\$^b\$")
                end
            else
                print(" & --\$^c\$")
            end
        end
        println(" \\\\")
    end

    # Add spacing between warehouse groups except for the last one
    if w_idx < length(warehouse_levels)
        println("\\midrule")
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println(
    "      \\item \\textit{Notes.} Split ratios are percentages (lower is better). Diff and Buffer columns show the size difference and buffer parameters respectively.",
)
println("      \\item \$^a\$ Reduced sample size due to timeouts on larger instances.")
println("      \\item \$^b\$ No solutions found within the time limit.")
println("      \\item \$^c\$ Algorithm not evaluated for this configuration.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")

# Generate detailed analysis
println("\n% Detailed Analysis by Configuration:")
println("% ===================================")

for wareh in warehouse_levels
    println("% Warehouse $wareh:")
    for (diff, buffer) in config_combinations
        config_results = results[
            (results.wareh .== wareh) .& (results.diff .== diff) .& (results.buffer .== buffer),
            :,
        ]

        if nrow(config_results) > 0
            println("%   Config ($(diff), $(buffer)):")
            for heur in heuristics
                heur_data = config_results[config_results.heuristic .== heur, :]
                if nrow(heur_data) > 0
                    ratio = heur_data[1, :split_ratio]
                    success_rate = heur_data[1, :success_rate]
                    if ratio !== nothing
                        println(
                            "%     $heur: $ratio% ($(round(success_rate, digits=1))% success)",
                        )
                    else
                        println(
                            "%     $heur: No solutions ($(round(success_rate, digits=1))% success)",
                        )
                    end
                else
                    println("%     $heur: Not tested")
                end
            end
        end
    end
    println("%")
end
