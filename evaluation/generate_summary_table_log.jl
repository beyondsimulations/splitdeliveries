using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Read the CSV file
df = CSV.read("results/overall_results.csv", DataFrame)

# Uniform-weight tables only; the weighted modes get dedicated tables (B3).
df = df[df.weight_mode .== "uniform", :]

# Success definition (B2): a run counts as successful if it returned a feasible
# allocation within the practical wall-clock cap. The solver time limit is
# 900 s; the cap adds tolerance for model building. DECISION PENDING
# (publication_checklist.md, B2): whether OPT should additionally require
# gap == 0 to be reported, and the exact cap value.
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
    "EMCI" => "MCI",
    "IIH" => "IIH",
    "IIHS" => "IIHS",
    "RND" => "RND",
)

println("=== DIAGNOSTIC ANALYSIS ===")
println("Total rows in dataset: $(nrow(df))")
println("Unique modes: $(sort(unique(df.mode)))")
println()

# Filter for relevant modes and map to standard names
filtered_df = df[in.(df.mode, Ref(keys(mode_mapping))), :]
filtered_df.heuristic = [mode_mapping[mode] for mode in filtered_df.mode]

println("After filtering and mapping:")
println("Rows: $(nrow(filtered_df))")
println("Heuristics: $(sort(unique(filtered_df.heuristic)))")
println()

# Add split ratio calculation
filtered_df.split_ratio = filtered_df.parcel_test ./ filtered_df.orders

# Check for any invalid split ratios
invalid_ratios = sum(isnan.(filtered_df.split_ratio) .| isinf.(filtered_df.split_ratio))
println("Invalid split ratios: $invalid_ratios")
println(
    "Split ratio range: $(minimum(filtered_df.split_ratio)) to $(maximum(filtered_df.split_ratio))",
)
println()

# Define heuristics in order for table
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "MCI", "IIH", "IIHS", "RND"]

println("=== SAMPLE SIZE ANALYSIS ===")
for heur in heuristics
    total_instances = nrow(filtered_df[filtered_df.heuristic .== heur, :])
    successful_instances = nrow(
        filtered_df[(filtered_df.heuristic .== heur) .& (filtered_df.duration .< SUCCESS_CAP), :]
    )
    timeout_instances = total_instances - successful_instances

    if total_instances > 0
        success_rate = (successful_instances / total_instances) * 100
        println(
            "$heur: $successful_instances/$total_instances successful ($(round(success_rate, digits=1))%)",
        )

        if timeout_instances > 0
            timeout_subset = filtered_df[
                (filtered_df.heuristic .== heur) .& (filtered_df.duration .>= SUCCESS_CAP), :,
            ]
            if nrow(timeout_subset) > 0
                println("  Timeouts by SKU: $(sort(unique(timeout_subset.skus)))")
            end
        end
    else
        println("$heur: No instances found")
    end
end
println()

println("=== SPLIT RATIO ANALYSIS ===")
# Calculate statistics for each heuristic
for heur in heuristics
    subset = filtered_df[
        (filtered_df.heuristic .== heur) .& (filtered_df.duration .< SUCCESS_CAP), :,
    ]

    if nrow(subset) > 0
        split_ratios = subset.split_ratio * 100  # Convert to percentage

        println("$heur (n=$(nrow(subset))):")
        println("  Mean: $(round(mean(split_ratios), digits=2))%")
        println("  Median: $(round(median(split_ratios), digits=2))%")
        println("  Min: $(round(minimum(split_ratios), digits=2))%")
        println("  Max: $(round(maximum(split_ratios), digits=2))%")
        println("  Std: $(round(std(split_ratios), digits=2))%")

        # Check distribution by dataset characteristics
        println("  By SKU levels:")
        for sku in sort(unique(subset.skus))
            sku_subset = subset[subset.skus .== sku, :]
            if nrow(sku_subset) > 0
                sku_mean = mean(sku_subset.split_ratio) * 100
                println(
                    "    SKU $sku: $(round(sku_mean, digits=2))% (n=$(nrow(sku_subset)))"
                )
            end
        end
        println()
    else
        println("$heur: No successful instances")
        println()
    end
end

println("=== INSTANCE-WISE IMPROVEMENT ANALYSIS ===")
# Get RND baseline
rnd_subset = filtered_df[
    (filtered_df.heuristic .== "RND") .& (filtered_df.duration .< SUCCESS_CAP), :,
]
if nrow(rnd_subset) > 0
    println("RND baseline: $(nrow(rnd_subset)) total instances")

    for heur in filter(h -> h != "RND", heuristics)
        subset = filtered_df[
            (filtered_df.heuristic .== heur) .& (filtered_df.duration .< SUCCESS_CAP), :,
        ]

        if nrow(subset) > 0
            # Instance-by-instance comparison
            # Create lookup for RND results by instance configuration
            rnd_lookup = Dict()
            for row in eachrow(rnd_subset)
                key = (row.skus, row.wareh, row.diff, row.buffer, row.dependency)
                rnd_lookup[key] = row.split_ratio
            end

            # Calculate instance-wise improvements
            improvements = Float64[]

            for row in eachrow(subset)
                key = (row.skus, row.wareh, row.diff, row.buffer, row.dependency)
                if haskey(rnd_lookup, key)
                    rnd_ratio = rnd_lookup[key]
                    heur_ratio = row.split_ratio
                    improvement = ((rnd_ratio - heur_ratio) / rnd_ratio) * 100
                    push!(improvements, improvement)
                end
            end

            if length(improvements) > 0
                println(
                    "$heur - Instance-wise improvements vs RND (n=$(length(improvements))):"
                )
                println("  Mean: $(round(mean(improvements), digits=2))%")
                println("  Median: $(round(median(improvements), digits=2))%")
                println("  Min: $(round(minimum(improvements), digits=2))%")
                println("  Max: $(round(maximum(improvements), digits=2))%")
                println("  Std: $(round(std(improvements), digits=2))%")
                println("  Common instances: $(length(improvements))/$(nrow(rnd_subset))")
            else
                println("$heur: No common instances with RND")
            end
        end
    end
else
    println("No RND instances found!")
end

println("\n=== BEST PERFORMER ANALYSIS ===")
# Find which heuristic achieves the best (lowest) split ratio for each instance
function analyze_best_performers(df)
    println("Analyzing which heuristic achieves the lowest split ratio per instance...")

    # Create DataFrame with all successful instances
    all_successful = df[df.duration .< SUCCESS_CAP, :]
    println("Total successful instance-heuristic combinations: $(nrow(all_successful))")

    # Group by instance configuration and find winner for each
    instance_groups = groupby(all_successful, [:skus, :wareh, :diff, :buffer, :dependency])
    winners = Dict{String,Float64}()

    # Count competing instances and winners
    num_competing = 0
    for group in instance_groups
        if nrow(group) > 1  # Only consider instances where multiple heuristics succeeded
            min_split_ratio = minimum(group.split_ratio)
            best_performers = group[group.split_ratio .== min_split_ratio, :]

            # In case of ties, count each tied heuristic
            for row in eachrow(best_performers)
                heuristic = row.heuristic
                credit = 1.0 / nrow(best_performers)  # Split credit for ties
                winners[heuristic] = get(winners, heuristic, 0.0) + credit
            end
            num_competing += 1
        end
    end

    println("Instances with multiple competing heuristics: $num_competing")
    println(
        "\nBest performer frequency (instances where this heuristic achieved lowest split ratio):",
    )

    # Sort heuristics by win count
    sorted_winners = sort(collect(winners); by = x->x[2], rev = true)

    for (heuristic, wins) in sorted_winners
        win_percentage = (wins / num_competing) * 100
        println(
            "  $heuristic: $(round(wins, digits=1))/$num_competing instances ($(round(win_percentage, digits=1))%)",
        )
    end

    println("\n=== INSTANCE PROPERTIES WHERE EACH HEURISTIC EXCELS ===")

    # Analyze properties of instances where each heuristic wins
    for group in instance_groups
        if nrow(group) > 1
            min_split_ratio = minimum(group.split_ratio)
            best_performers = group[group.split_ratio .== min_split_ratio, :]

            # Store winning instances for each heuristic
            for row in eachrow(best_performers)
                heuristic = row.heuristic
                if !haskey(winners, heuristic)
                    continue
                end

                # Create a key for this winning instance
                instance_key = (
                    heuristic = heuristic,
                    skus = row.skus,
                    wareh = row.wareh,
                    diff = row.diff,
                    buffer = row.buffer,
                    dependency = row.dependency,
                    split_ratio = row.split_ratio,
                )

                # We'll collect these for analysis
            end
        end
    end

    # Now analyze properties for top performers
    top_heuristics = [h for (h, _) in sorted_winners if get(winners, h, 0) >= 10.0]  # At least 10 wins

    for heuristic in top_heuristics
        println("\n$heuristic - Properties of best-performing instances:")

        # Collect all instances where this heuristic wins
        winning_instances = []
        for group in instance_groups
            if nrow(group) > 1
                min_split_ratio = minimum(group.split_ratio)
                best_performers = group[group.split_ratio .== min_split_ratio, :]

                for row in eachrow(best_performers)
                    if row.heuristic == heuristic
                        push!(winning_instances, row)
                    end
                end
            end
        end

        if length(winning_instances) > 0
            println("  SKU distribution:")
            sku_counts = Dict()
            for instance in winning_instances
                sku_counts[instance.skus] = get(sku_counts, instance.skus, 0) + 1
            end
            for (sku, count) in sort(collect(sku_counts))
                pct = (count / length(winning_instances)) * 100
                println("    SKU $sku: $count instances ($(round(pct, digits=1))%)")
            end

            println("  Warehouse distribution:")
            wareh_counts = Dict()
            for instance in winning_instances
                wareh_counts[instance.wareh] = get(wareh_counts, instance.wareh, 0) + 1
            end
            for (wareh, count) in sort(collect(wareh_counts))
                pct = (count / length(winning_instances)) * 100
                println(
                    "    $wareh warehouses: $count instances ($(round(pct, digits=1))%)"
                )
            end

            println("  Difficulty distribution:")
            diff_counts = Dict()
            for instance in winning_instances
                diff_counts[instance.diff] = get(diff_counts, instance.diff, 0) + 1
            end
            for (diff, count) in sort(collect(diff_counts))
                pct = (count / length(winning_instances)) * 100
                println("    Difficulty $diff: $count instances ($(round(pct, digits=1))%)")
            end

            println("  Buffer distribution:")
            buffer_counts = Dict()
            for instance in winning_instances
                buffer_counts[instance.buffer] = get(buffer_counts, instance.buffer, 0) + 1
            end
            for (buffer, count) in sort(collect(buffer_counts))
                pct = (count / length(winning_instances)) * 100
                println("    Buffer $buffer: $count instances ($(round(pct, digits=1))%)")
            end

            println("  Dependency distribution:")
            dep_counts = Dict()
            for instance in winning_instances
                dep_counts[instance.dependency] =
                    get(dep_counts, instance.dependency, 0) + 1
            end
            for (dep, count) in sort(collect(dep_counts))
                pct = (count / length(winning_instances)) * 100
                println("    Dependency $dep: $count instances ($(round(pct, digits=1))%)")
            end

            # Show best and worst performance
            split_ratios = [instance.split_ratio * 100 for instance in winning_instances]
            println("  Performance on winning instances:")
            println("    Best split ratio: $(round(minimum(split_ratios), digits=2))%")
            println("    Worst split ratio: $(round(maximum(split_ratios), digits=2))%")
            println("    Mean split ratio: $(round(mean(split_ratios), digits=2))%")
        end
    end
end

analyze_best_performers(filtered_df)

println("\n=== INSTANCE OVERLAP ANALYSIS ===")
# Check what instances each heuristic successfully solves
instance_keys = []
for row in eachrow(filtered_df)
    if row.duration < SUCCESS_CAP
        push!(
            instance_keys,
            (
                heuristic = row.heuristic,
                skus = row.skus,
                wareh = row.wareh,
                diff = row.diff,
                buffer = row.buffer,
                dependency = row.dependency,
            ),
        )
    end
end

instance_df = DataFrame(instance_keys)
println("Total successful instance-heuristic combinations: $(nrow(instance_df))")

# Group by instance configuration
instance_configs = groupby(instance_df, [:skus, :wareh, :diff, :buffer, :dependency])
println("Unique instance configurations: $(length(instance_configs))")

# Check how many heuristics solve each configuration
config_coverage = []
for config in instance_configs
    heuristics_for_config = unique(config.heuristic)
    push!(config_coverage, length(heuristics_for_config))
end

println(
    "Heuristics per configuration - Min: $(minimum(config_coverage)), Max: $(maximum(config_coverage)), Mean: $(round(mean(config_coverage), digits=1))",
)
