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

# Add split ratio calculation for test
filtered_df.split_ratio = filtered_df.parcel_test ./ filtered_df.orders

# Define heuristics in order for table
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "RND"]

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

# Calculate results for each dependency structure and heuristic
results = DataFrame(;
    dependency = String[], heuristic = String[], split_ratio = Any[], sample_size = Int[]
)

println("Dataset structure analysis:")
println("="^60)

for dep in dependency_levels
    println("\nDependency structure: $dep")

    for heur in heuristics
        subset = filtered_df[
            (filtered_df.dependency .== dep) .& (filtered_df.heuristic .== heur), :,
        ]

        if nrow(subset) > 0
            # Calculate success rate and average split ratio
            successful_runs = subset.duration[subset.duration .< SUCCESS_CAP]
            successful_ratios = subset.split_ratio[subset.duration .< SUCCESS_CAP]
            total_runs = nrow(subset)
            success_rate = length(successful_runs) / total_runs * 100

            println(
                "  $heur: $(length(successful_ratios))/$(total_runs) successful ($(round(success_rate, digits=1))%)",
            )

            if length(successful_ratios) > 0
                # Calculate actual split ratio as percentage
                avg_ratio = mean(successful_ratios)
                split_ratio = round(avg_ratio * 100; digits = 2)
                sample_size = length(successful_ratios)
            else
                split_ratio = nothing
                sample_size = 0
            end
        else
            split_ratio = nothing
            sample_size = 0
        end

        # Add to results if we have a valid result
        if split_ratio !== nothing
            push!(
                results,
                (
                    dependency = dep,
                    heuristic = heur,
                    split_ratio = split_ratio,
                    sample_size = sample_size,
                ),
            )
        end
    end
end

# Display the table
println("\n\n")
println("\\caption{Average split ratio (in \\%) depending on the dataset structure}")
println("\\label{tab:dep}")
println("\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrr}")
println("\\toprule")
print("\t \t")
for heur in heuristics
    print(" & $heur")
end
println("\\\\")
println("\\midrule")

# Print data rows
for (i, dep) in enumerate(dependency_levels)
    dep_results = results[results.dependency .== dep, :]

    # Format dependency name
    dep_formatted = dep * "\$^b\$"

    # Split ratios
    print("$dep_formatted")
    for heur in heuristics
        heur_data = dep_results[dep_results.heuristic .== heur, :]
        if nrow(heur_data) > 0
            ratio = heur_data[1, :split_ratio]
            sample_size = heur_data[1, :sample_size]

            # Check if this algorithm has reduced sample size compared to others
            all_samples = [row.sample_size for row in eachrow(dep_results)]
            max_sample = maximum(all_samples)
            has_reduced_sample = sample_size < max_sample * 0.9  # Less than 90% of max

            if has_reduced_sample && heur != "RND"
                print(" & $ratio\$^a\$")
            else
                print(" & $ratio")
            end
        else
            print(" & -")
        end
    end
    println(" \\\\")

    # Add midrule after ID and MD groups
    if i == 2 || i == 4
        println("\\midrule")
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println(
    "      \\item \\textit{Notes.} The split ratio is calculated using the number of split parcels and the number of orders.",
)
println(
    "      \\item \$^a\$ Smaller sample size, as some larger instances could not be solved."
)
println("      \\item \$^b\$ Further details are described in Table \\ref{tab:datasets}.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")
