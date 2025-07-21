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
println("Split ratio range: $(minimum(filtered_df.split_ratio)) to $(maximum(filtered_df.split_ratio))")
println()

# Define heuristics in order for table
heuristics = ["QMK", "CHI", "KL", "GO", "GP", "GS", "BS", "RND"]

println("=== SAMPLE SIZE ANALYSIS ===")
for heur in heuristics
    total_instances = nrow(filtered_df[filtered_df.heuristic.==heur, :])
    successful_instances = nrow(filtered_df[(filtered_df.heuristic.==heur).&(filtered_df.duration.<3600), :])
    timeout_instances = total_instances - successful_instances

    if total_instances > 0
        success_rate = (successful_instances / total_instances) * 100
        println("$heur: $successful_instances/$total_instances successful ($(round(success_rate, digits=1))%)")

        if timeout_instances > 0
            timeout_subset = filtered_df[(filtered_df.heuristic.==heur).&(filtered_df.duration.>=3600), :]
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
    subset = filtered_df[(filtered_df.heuristic.==heur).&(filtered_df.duration.<3600), :]

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
            sku_subset = subset[subset.skus.==sku, :]
            if nrow(sku_subset) > 0
                sku_mean = mean(sku_subset.split_ratio) * 100
                println("    SKU $sku: $(round(sku_mean, digits=2))% (n=$(nrow(sku_subset)))")
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
rnd_subset = filtered_df[(filtered_df.heuristic.=="RND").&(filtered_df.duration.<3600), :]
if nrow(rnd_subset) > 0
    println("RND baseline: $(nrow(rnd_subset)) total instances")

    for heur in filter(h -> h != "RND", heuristics)
        subset = filtered_df[(filtered_df.heuristic.==heur).&(filtered_df.duration.<3600), :]

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
                println("$heur - Instance-wise improvements vs RND (n=$(length(improvements))):")
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

println("\n=== INSTANCE OVERLAP ANALYSIS ===")
# Check what instances each heuristic successfully solves
instance_keys = []
for row in eachrow(filtered_df)
    if row.duration < 3600
        push!(instance_keys, (heuristic=row.heuristic,
            skus=row.skus, wareh=row.wareh, diff=row.diff,
            buffer=row.buffer, dependency=row.dependency))
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

println("Heuristics per configuration - Min: $(minimum(config_coverage)), Max: $(maximum(config_coverage)), Mean: $(round(mean(config_coverage), digits=1))")
