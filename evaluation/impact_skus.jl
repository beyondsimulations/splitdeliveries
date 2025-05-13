using DataFrames, CSV, Statistics
using CairoMakie  # For publication-quality plots
using ColorSchemes  # For better color schemes

# Read the overall results
df = CSV.read("results/overall_results.csv", DataFrame)

# Filter for the modes we're interested in
modes_of_interest = ["CHI", "QMK", "KL", "GO", "GP", "GS", "BS"]

# Create summary statistics grouped by SKUs, orders, and mode
# Filter out failed runs (where duration is 0 or parcel_test is 0)
summary_stats = combine(groupby(filter(row -> row.duration > 0 && row.parcel_test > 0, df), [:skus, :orders, :mode]), 
    :duration => mean => :avg_duration,
    :parcel_test => mean => :avg_parcels,
    :parcel_train => mean => :avg_parcels_train)

# Calculate improvement over random allocation
function calculate_improvement(df, metric)
    improvements = DataFrame()
    
    # Get random baseline for each configuration
    random_data = filter(row -> row.mode == "RND" && row.duration > 0 && row.parcel_test > 0, df)
    random_means = combine(groupby(random_data, [:skus, :orders]), 
        Symbol(metric) => mean => :random_mean)
    
    # Calculate improvement for each algorithm
    for mode in modes_of_interest
        mode_data = filter(row -> row.mode == mode && row.duration > 0 && row.parcel_test > 0, df)
        mode_means = combine(groupby(mode_data, [:skus, :orders]), 
            Symbol(metric) => mean => :algo_mean)
        
        # Join with random baseline
        joined = leftjoin(mode_means, random_means, on=[:skus, :orders])
        
        # Calculate improvement percentage
        joined.improvement = ((joined.random_mean .- joined.algo_mean) ./ joined.random_mean) .* 100
        
        # Add mode column correctly
        joined[!, :mode] .= mode
        
        append!(improvements, joined)
    end
    return improvements
end

# Calculate improvements for parcels
parcel_improvements = calculate_improvement(df, :parcel_test)

# Print success rates
println("\nAlgorithm Success Rates:")
println("----------------------")
for mode in modes_of_interest
    total_runs = nrow(filter(row -> row.mode == mode, df))
    successful_runs = nrow(filter(row -> row.mode == mode && row.duration > 0 && row.parcel_test > 0, df))
    success_rate = (successful_runs / total_runs) * 100
    println("$mode: $(round(success_rate, digits=2))% ($successful_runs/$total_runs)")
end

# Set up the color scheme
colors = ColorSchemes.tab10.colors
mode_colors = Dict(mode => colors[i] for (i, mode) in enumerate(modes_of_interest))

# Create a figure with subplots
fig = Figure(size=(1200, 1000))

# 1. Improvement in Split Deliveries vs Number of Orders
ax1 = Axis(fig[1, 1], 
    xlabel="Number of Orders", 
    ylabel="Improvement in Split Deliveries (%)",
    title="Improvement over Random Allocation by Algorithm")

for mode in modes_of_interest
    mode_data = filter(row -> row.mode == mode, parcel_improvements)
    if !isempty(mode_data)
        scatter!(ax1, mode_data.orders, mode_data.improvement, 
            color=mode_colors[mode],
            label=mode,
            markersize=10)
    end
end

# 2. Computational Time vs Number of Orders
ax2 = Axis(fig[1, 2], 
    xlabel="Number of Orders", 
    ylabel="Average Computational Time (seconds)",
    title="Computational Time by Algorithm")

for mode in modes_of_interest
    mode_data = filter(row -> row.mode == mode, summary_stats)
    if !isempty(mode_data)
        scatter!(ax2, mode_data.orders, mode_data.avg_duration, 
            color=mode_colors[mode],
            label=mode,
            markersize=10)
    end
end

# 3. Improvement in Split Deliveries by Dependency Type
ax3 = Axis(fig[2, 1], 
    xlabel="Dependency Type", 
    ylabel="Improvement in Split Deliveries (%)",
    title="Improvement over Random by Dependency Type")

# Calculate average improvement by dependency type
dep_improvements = combine(groupby(parcel_improvements, [:mode]), 
    :improvement => mean => :avg_improvement)

for (i, mode) in enumerate(modes_of_interest)
    mode_data = filter(row -> row.mode == mode, dep_improvements)
    if !isempty(mode_data)
        barplot!(ax3, [i], [mode_data.avg_improvement[1]], 
            color=mode_colors[mode],
            label=mode)
    end
end

# 4. Improvement vs Computational Time Trade-off
ax4 = Axis(fig[2, 2], 
    xlabel="Computational Time (seconds)", 
    ylabel="Improvement in Split Deliveries (%)",
    title="Improvement vs Computational Time Trade-off")

for mode in modes_of_interest
    mode_data = filter(row -> row.mode == mode, summary_stats)
    imp_data = filter(row -> row.mode == mode, parcel_improvements)
    if !isempty(mode_data) && !isempty(imp_data)
        scatter!(ax4, mode_data.avg_duration, imp_data.improvement, 
            color=mode_colors[mode],
            label=mode,
            markersize=10)
    end
end

# Add a single legend for all plots
Legend(fig[3, :], ax1, "Algorithms", orientation=:horizontal)

# Adjust layout
colgap!(fig.layout, 20)
rowgap!(fig.layout, 20)

# Save the figure
save("results/algorithm_improvements.png", fig, px_per_unit=2)

# Print detailed summary statistics
println("\nDetailed Summary Statistics by Algorithm:")
println("----------------------------------------")
for mode in modes_of_interest
    mode_data = filter(row -> row.mode == mode, parcel_improvements)
    if !isempty(mode_data)
        println("\nAlgorithm: $mode")
        println("Average Improvement over Random: $(round(mean(mode_data.improvement), digits=2))%")
        println("Standard Deviation of Improvement: $(round(std(mode_data.improvement), digits=2))%")
        
        # Print statistics by order count
        println("\nBy Order Count:")
        for order in unique(mode_data.orders)
            order_data = filter(row -> row.orders == order, mode_data)
            println("  Orders: $order")
            println("    Average Improvement: $(round(mean(order_data.improvement), digits=2))%")
        end
    end
end

# Print computational time statistics
println("\nComputational Time Statistics:")
println("----------------------------")
for mode in modes_of_interest
    mode_data = filter(row -> row.mode == mode, summary_stats)
    if !isempty(mode_data)
        println("\nAlgorithm: $mode")
        println("Average Computational Time: $(round(mean(mode_data.avg_duration), digits=2)) seconds")
        println("Standard Deviation: $(round(std(mode_data.avg_duration), digits=2)) seconds")
    end
end