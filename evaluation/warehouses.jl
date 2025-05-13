using DataFrames, CSV, Statistics, Printf

# Read the overall results
df = CSV.read("results/overall_results.csv", DataFrame)

# Define the algorithms we want to include
algorithms = ["OPT", "QMK", "QMKJ", "CHI", "CHIM", "KL", "GO", "GP", "GS", "BS", "RND"]

# Get all unique warehouse-buffer-diff combinations from the data
df_combinations = unique(df[:, [:wareh, :buffer, :diff]])
sort!(df_combinations, [:wareh, :buffer, :diff])

# Get total number of instances for each combination from RND mode
total_instances = combine(
    groupby(filter(row -> row.mode == "RND", df), [:wareh, :buffer, :diff]), 
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

# Create DataFrames for calculations
split_ratio_by_combo = DataFrame(wareh = df_combinations.wareh, 
                                buffer = df_combinations.buffer,
                                diff = df_combinations.diff)

# Calculate split ratios for each algorithm
for algo in algorithms
    split_ratio_by_combo[!, algo] = [
        begin
            algo_data = filter(
                row -> row.mode == algo && 
                       row.wareh == wareh && 
                       row.buffer == buffer &&
                       row.diff == diff && 
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
        for (wareh, buffer, diff) in zip(df_combinations.wareh, df_combinations.buffer, df_combinations.diff)
    ]
end

# Format the results for LaTeX
println("\n\\begin{threeparttable}")
println("\\begin{tabular}{lllrrrrrrrrrrr}")
println("\\toprule")
println("\$|W|\$ & \$B\$ & \$D\$ & OPT & QMK & QMKJ & CHI & CHIM & KL & GO & GP & GS & BS & RND \\\\")
println("\\midrule")

# Declare current_wareh as a local variable to avoid scope issues
let current_wareh = -1
    for i in 1:nrow(df_combinations)
        wareh = df_combinations.wareh[i]
        buffer = df_combinations.buffer[i]
        diff = df_combinations.diff[i]
        
        # Add extra spacing between different warehouse groups
        if wareh != current_wareh && i > 1
            println("\\midrule")
            current_wareh = wareh
        elseif i == 1
            current_wareh = wareh
        end
        
        wareh_str = @sprintf("%d", Int(wareh))
        buffer_str = @sprintf("%.1f", buffer)
        diff_str = @sprintf("%.1f", diff)
        
        # Format split ratio values and success rates
        ratios = []
        rates = []
        
        for algo in algorithms
            # Find matching total instances
            total_row = filter(r -> r.wareh == wareh && r.buffer == buffer && r.diff == diff, total_instances)
            if isempty(total_row)
                total = 0
            else
                total = total_row[1, :total]
            end
            
            # Get success rate
            successful = nrow(filter(
                r -> r.mode == algo && 
                     r.wareh == wareh && 
                     r.buffer == buffer &&
                     r.diff == diff && 
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
        ratio_row = join([wareh_str, buffer_str, diff_str, ratios...], " & ")
        rate_row = join(["", "", "", rates...], " & ")
        println("$ratio_row \\\\")
        println("$rate_row \\\\")
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println("      \\item \\textit{Notes.} The split ratio is displayed in the first row, with the success rate shown in the second row. Split ratio is calculated as the number of parcels divided by (number of orders × train_test fraction). \$|W|\$ represents the number of warehouses, \$B\$ represents the buffer, and \$D\$ represents the difference. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.")
println("      \\item \$^a\$ No solution could be found within 3600 seconds.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")

# Print raw data for verification
println("\nRaw split ratio data:")
println(split_ratio_by_combo)
