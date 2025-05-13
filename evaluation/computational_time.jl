using DataFrames, CSV, Statistics, Printf

# Read the overall results
df = CSV.read("results/overall_results.csv", DataFrame)

# Define the algorithms we want to include
algorithms = ["OPT","QMK", "QMKJ", "CHI", "CHIM", "KL", "GO", "GP", "GS", "BS"]

# Get all unique SKU numbers from the data
# Get all unique SKU numbers from the data
sku_ranges = sort(unique(df.skus))

# First, let's analyze the successful instances (time <= 900 seconds)
println("Success Analysis (instances that completed within 900 seconds):")
println("--------------------------------------------------")

# Get total number of instances for each SKU count from RND mode
total_instances = combine(groupby(filter(row -> row.mode == "RND", df), :skus), nrow => :total)

# Analyze successful instances for each algorithm and SKU count
for algo in algorithms
    println("\nAlgorithm: $algo")
    println("SKU count | Successful | Total instances | Success rate")
    println("------------------------------------------------")
    
    for sku in sku_ranges
        # Get total instances for this SKU count
        total = total_instances[total_instances.skus .== sku, :total][1]
        
        # Count successful instances (duration <= 900)
        successful = nrow(filter(row -> row.mode == algo && row.skus == sku && row.duration <= 920, df))
        
        # Calculate success rate
        success_rate = (successful / total) * 100
        
        println("$sku | $successful | $total | $(round(Int, success_rate))%")
    end
end

# Function to format time values
function format_time(time, success_rate)
    if time >= 920 || (time == 0.0 && success_rate == 0)  # If time is greater than 1 hour or time is 0 with 0% success
        return "\$^a\$"  # Use the footnote marker
    else
        return @sprintf("%.2f", round(time, digits=2))  # Round to 2 decimal places
    end
end

# Create a DataFrame to store the results
results = DataFrame(
    skus = sku_ranges,
    OPT = zeros(length(sku_ranges)),
    QMK = zeros(length(sku_ranges)),
    QMKJ = zeros(length(sku_ranges)),
    CHI = zeros(length(sku_ranges)),
    CHIM = zeros(length(sku_ranges)),
    KL = zeros(length(sku_ranges)),
    GO = zeros(length(sku_ranges)),
    GP = zeros(length(sku_ranges)),
    GS = zeros(length(sku_ranges)),
    BS = zeros(length(sku_ranges))
)

# Calculate average times for each SKU range and algorithm
for (i, sku) in enumerate(sku_ranges)
    for algo in algorithms
        # Filter data for this SKU range and algorithm
        algo_data = filter(row -> row.mode == algo && row.skus == sku && row.duration > 0 && row.duration < 950, df)
        
        if !isempty(algo_data)
            # Calculate mean duration and round to 2 decimal places
            avg_time = round(mean(algo_data.duration), digits=2)
            results[i, Symbol(algo)] = avg_time
        end
    end
end

# Format the results for LaTeX
println("\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrr}")
println("\\toprule")
println("\$|\\mathcal{I}|\$  & OPT & QMK & QMKJ & CHI & CHIM & KL & GO & GP & GS & BS \\\\")
println("\\midrule")

for row in eachrow(results)
    # Format SKU count with comma
    sku_str = @sprintf("%d", Int(row.skus))
    
    # Format computation times and success rates
    times = []
    rates = []
    for algo in algorithms
        # Get total instances for this SKU count
        total = total_instances[total_instances.skus .== row.skus, :total][1]
        # Count successful instances for this specific algorithm and SKU count
        successful = nrow(filter(r -> r.mode == algo && r.skus == row.skus && r.duration <= 920, df))
        success_rate = round(Int, (successful / total) * 100)
        
        # Format time based on time and success rate
        time_str = format_time(row[Symbol(algo)], success_rate)
        
        push!(times, "\$$time_str\$")
        
        # Skip showing 0% success rates when using the footnote marker
        if time_str == "\$^a\$"
            push!(rates, "")
        else
            push!(rates, "\$$success_rate\\%\$")
        end
    end
    
    # Create the rows for times and rates
    time_row = join([sku_str, times...], " & ")
    rate_row = join(["", rates...], " & ")
    println("$time_row \\\\")
    println("$rate_row \\\\")
    println("\\midrule")
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println("      \\item \\textit{Notes.} The computation time is displayed in seconds in the first row, with the success rate shown in the second row. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.")
println("      \\item \$^a\$ No solution could be found within 900 seconds.")

println("\\end{tablenotes}")
println("\\end{threeparttable}")

# Also print the raw data for verification
println("\nRaw data for verification:")
println(results)
