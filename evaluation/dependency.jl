using DataFrames, CSV, Statistics, Printf

# Read the overall results
df = CSV.read("results/overall_results.csv", DataFrame)

# Define the algorithms we want to include
algorithms = ["OPT", "QMK", "QMKJ", "CHI", "CHIM", "KL", "GO", "GP", "GS", "BS", "RND"]
algorithms_no_rnd = filter(a -> a != "RND", algorithms)

# Get all unique dependency combinations from the data
df_combinations = unique(df[:, [:dependency]])
sort!(df_combinations, [:dependency])

# Get total number of instances for each dependency from RND mode
total_instances = combine(
    groupby(filter(row -> row.mode == "RND", df), [:dependency]), 
    nrow => :total
)

# Function to format split ratio values
function format_split_ratio(ratio, success_rate)
    if isnan(ratio) || success_rate == 0  # If ratio is NaN or no successful instances
        return "\$^a\$"  # Use the footnote marker
    else
        return @sprintf("%.2f", ratio)  # Round to 2 decimal places
    end
end

# Function to format percentage reduction
function format_reduction(reduction, success_rate)
    if isnan(reduction) || success_rate == 0  # If reduction is NaN or no successful instances
        return "\$^a\$"  # Use the footnote marker
    else
        return @sprintf("%.1f\\%%", reduction * 100)  # Format as percentage
    end
end

# Create DataFrames for calculations
split_ratio_by_combo = DataFrame(dependency = df_combinations.dependency)

# Calculate split ratios for each algorithm
for algo in algorithms
    split_ratio_by_combo[!, algo] = [
        begin
            algo_data = filter(
                row -> row.mode == algo && 
                       row.dependency == dependency && 
                       row.duration > 0 && 
                       row.duration < 3600, 
                df
            )
            if isempty(algo_data)
                NaN
            else
                # Calculate split ratio: parcels / (orders * train_test)
                mean(algo_data.parcel_test ./ (algo_data.orders .* algo_data.train_test))
            end
        end
        for dependency in df_combinations.dependency
    ]
end

# Create DataFrame for reductions compared to random
reduction_by_combo = DataFrame(dependency = df_combinations.dependency)

# Calculate reductions for each algorithm compared to RND
for algo in algorithms_no_rnd
    reduction_by_combo[!, algo] = [
        begin
            rnd_value = split_ratio_by_combo[i, "RND"]
            algo_value = split_ratio_by_combo[i, algo]
            if isnan(rnd_value) || isnan(algo_value) || rnd_value == 0
                NaN
            else
                # Calculate reduction: (RND - algo) / RND
                (rnd_value - algo_value) / rnd_value
            end
        end
        for i in 1:nrow(split_ratio_by_combo)
    ]
end

# Format the results for LaTeX - First table (split ratios)
println("\n\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrrr}")
println("\\toprule")
println("Dependency & OPT & QMK & QMKJ & CHI & CHIM & KL & GO & GP & GS & BS & RND \\\\")
println("\\midrule")

for i in 1:nrow(df_combinations)
    dependency = df_combinations.dependency[i]
    
    # Format split ratio values and success rates
    ratios = []
    rates = []
    
    for algo in algorithms
        # Find matching total instances
        total_row = filter(r -> r.dependency == dependency, total_instances)
        if isempty(total_row)
            total = 0
        else
            total = total_row[1, :total]
        end
        
        # Get success rate
        successful = nrow(filter(
            r -> r.mode == algo && 
                 r.dependency == dependency && 
                 r.duration <= 920, 
            df
        ))
        
        success_rate = total > 0 ? round(Int, (successful / total) * 100) : 0
        
        # Format split ratio
        ratio_value = split_ratio_by_combo[i, algo]
        ratio_str = format_split_ratio(ratio_value, success_rate)
        
        push!(ratios, "\$$ratio_str\$")
        
        # Skip showing 0% success rates when using the footnote marker
        if ratio_str == "\$^a\$"
            push!(rates, "")
        else
            push!(rates, "\$$success_rate\\%\$")
        end
    end
    
    # Create the rows for split ratios and rates
    ratio_row = join([dependency, ratios...], " & ")
    rate_row = join(["", rates...], " & ")
    println("$ratio_row \\\\")
    println("$rate_row \\\\")
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println("      \\item \\textit{Notes.} The split ratio is displayed in the first row, with the success rate shown in the second row. Split ratio is calculated as the number of parcels divided by (number of orders × train_test fraction). Dependency types: HD = High Dependency, MD = Medium Dependency, ID = Independent, SF = Same Frequency, VF = Variable Frequency. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.")
println("      \\item \$^a\$ No solution could be found within 3600 seconds.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")

# Format the results for LaTeX - Second table (reductions compared to random)
println("\n\\begin{threeparttable}")
println("\\begin{tabular}{lrrrrrrrrrr}")
println("\\toprule")
println("Dependency & OPT & QMK & QMKJ & CHI & CHIM & KL & GO & GP & GS & BS \\\\")
println("\\midrule")

for i in 1:nrow(df_combinations)
    dependency = df_combinations.dependency[i]
    
    # Format reduction values
    reductions = []
    
    for algo in algorithms_no_rnd
        # Find matching total instances
        total_row = filter(r -> r.dependency == dependency, total_instances)
        if isempty(total_row)
            total = 0
        else
            total = total_row[1, :total]
        end
        
        # Get success rate
        successful = nrow(filter(
            r -> r.mode == algo && 
                 r.dependency == dependency && 
                 r.duration <= 920, 
            df
        ))
        
        success_rate = total > 0 ? round(Int, (successful / total) * 100) : 0
        
        # Format reduction
        reduction_value = reduction_by_combo[i, algo]
        reduction_str = format_reduction(reduction_value, success_rate)
        
        push!(reductions, "\$$reduction_str\$")
    end
    
    # Create the row for reductions
    reduction_row = join([dependency, reductions...], " & ")
    println("$reduction_row \\\\")
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println("      \\item \\textit{Notes.} This table shows the percentage reduction in split ratio compared to random allocation. Higher percentages indicate better performance. Calculated as (RND - algo)/RND × 100%.")
println("      \\item \$^a\$ No reduction could be calculated due to missing data or no successful runs.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")

# Print raw data for verification
println("\nRaw split ratio data:")
println(split_ratio_by_combo)
println("\nRaw reduction data:")
println(reduction_by_combo)
