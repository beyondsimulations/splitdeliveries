using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Read the CSV file
df = CSV.read("results/overall_results.csv", DataFrame)

# Weighted-mode tables (frequency-proportional and random storage
# requirements); the uniform mode is covered by the dedicated tables.
weight_modes = ["frequency", "random"]
weight_labels = Dict(
    "frequency" => "Frequency-proportional storage requirements",
    "random" => "Random storage requirements",
)

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
    "QMKJ" => "QMK",
    "CHI_1.0e-5" => "CHI",
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

# Add split ratio calculation for test
filtered_df.split_ratio = filtered_df.parcel_test ./ filtered_df.orders

# Define heuristics in order for table
heuristics = ["QMK", "CHI", "GO", "GP", "GS", "BS", "EMCI", "RND"]

# Get unique dependency structures in desired order (ID, MD, HD)
all_deps = unique(filtered_df.dependency)
dependency_levels = []
for prefix in ["ID", "MD", "HD"]
    for dep in sort(all_deps)
        if startswith(dep, prefix)
            push!(dependency_levels, dep)
        end
    end
end

# Calculate results for each weight mode, dependency structure, and heuristic
results = DataFrame(;
    weight_mode = String[],
    dependency = String[],
    heuristic = String[],
    split_ratio = Any[],
    sample_size = Int[],
)

println("Weighted-mode analysis:")
println("="^60)

for wm in weight_modes
    println("\nWeight mode: $wm")
    for dep in dependency_levels
        for heur in heuristics
            subset = filtered_df[
                (filtered_df.weight_mode .== wm) .& (filtered_df.dependency .== dep) .& (filtered_df.heuristic .== heur),
                :,
            ]

            if nrow(subset) > 0
                success_mask = is_success.(eachrow(subset))
                successful_ratios = subset.split_ratio[success_mask]
                if length(successful_ratios) > 0
                    push!(
                        results,
                        (
                            weight_mode = wm,
                            dependency = dep,
                            heuristic = heur,
                            split_ratio = round(mean(successful_ratios) * 100; digits = 2),
                            sample_size = length(successful_ratios),
                        ),
                    )
                end
            end
        end
    end
end

# Display the table
println("\n\n")
println("\\caption{Average split ratio (in \\%) under heterogeneous storage requirements}")
println("\\label{tab:weighted}")
println("\\begin{threeparttable}")
println("\\begin{tabular}{l" * "r"^length(heuristics) * "}")
println("\\toprule")
print("\t \t")
for heur in heuristics
    print(" & $heur")
end
println("\\\\")

for wm in weight_modes
    println("\\midrule")
    println("\\multicolumn{$(length(heuristics) + 1)}{l}{\\textit{$(weight_labels[wm])}} \\\\")
    println("\\midrule")

    for (i, dep) in enumerate(dependency_levels)
        dep_results = results[
            (results.weight_mode .== wm) .& (results.dependency .== dep), :,
        ]
        print("$dep\$^b\$")
        max_sample =
            nrow(dep_results) > 0 ? maximum(dep_results.sample_size) : 0
        for heur in heuristics
            heur_data = dep_results[dep_results.heuristic .== heur, :]
            if nrow(heur_data) > 0
                ratio = heur_data[1, :split_ratio]
                # Mark columns evaluated on fewer instances (QMK scales only
                # to 1,000 SKUs)
                if heur_data[1, :sample_size] < max_sample * 0.9 && heur != "RND"
                    print(" & $ratio\$^a\$")
                else
                    print(" & $ratio")
                end
            else
                print(" & -")
            end
        end
        println(" \\\\")
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println(
    "      \\item \\textit{Notes.} The split ratio is calculated using the number of split parcels and the number of orders. Both modes run on buffered configurations, which KL cannot exploit; it is therefore not part of this comparison.",
)
println(
    "      \\item \$^a\$ Smaller sample size, as some larger instances could not be solved."
)
println("      \\item \$^b\$ Further details are described in Table \\ref{tab:datasets}.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")
