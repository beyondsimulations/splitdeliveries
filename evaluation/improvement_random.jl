using DataFrames, CSV, Statistics, Printf

# Read the overall results
df = CSV.read("results/overall_results.csv", DataFrame)

# Define the algorithms we want to include
algorithms = ["OPT", "QMK", "QMKJ", "CHI", "CHIM", "KL", "GO", "GP", "GS", "BS", "RND"]

# Get all unique SKU-order combinations from the data
df_combinations = unique(df[:, [:skus, :orders]])
sort!(df_combinations, [:skus, :orders])

# Get total number of instances for each SKU-order combination from RND mode
total_instances = combine(
    groupby(filter(row -> row.mode == "RND", df), [:skus, :orders]), 
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
split_ratio_by_combo = DataFrame(skus = df_combinations.skus, orders = df_combinations.orders)

# Calculate split ratios for each algorithm
for algo in algorithms
    split_ratio_by_combo[!, algo] = [
        begin
            algo_data = filter(
                row -> row.mode == algo && 
                       row.skus == sku && 
                       row.orders == order_size && 
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
        for (sku, order_size) in zip(df_combinations.skus, df_combinations.orders)
    ]
end

# Format the results for LaTeX
println("\n\\begin{threeparttable}")
println("\\begin{tabular}{llrrrrrrrrrrr}")
println("\\toprule")
println("\$|\\mathcal{I}|\$ & \$|O|\$ & OPT & QMK & QMKJ & CHI & CHIM & KL & GO & GP & GS & BS & RND \\\\")
println("\\midrule")

# Declare current_sku as a local variable to avoid scope issues
let current_sku = -1
    for i in 1:nrow(df_combinations)
        sku = df_combinations.skus[i]
        order_size = df_combinations.orders[i]
        
        # Add extra spacing between different SKU groups
        if sku != current_sku && i > 1
            println("\\midrule")
            current_sku = sku
        elseif i == 1
            current_sku = sku
        end
        
        sku_str = @sprintf("%d", Int(sku))
        order_str = @sprintf("%d", Int(order_size/sku))
        
        # Format split ratio values and success rates
        ratios = []
        rates = []
        
        for algo in algorithms
            # Find matching total instances
            total_row = filter(r -> r.skus == sku && r.orders == order_size, total_instances)
            if isempty(total_row)
                total = 0
            else
                total = total_row[1, :total]
            end
            
            # Get success rate
            successful = nrow(filter(
                r -> r.mode == algo && 
                     r.skus == sku && 
                     r.orders == order_size && 
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
        ratio_row = join([sku_str, order_str, ratios...], " & ")
        rate_row = join(["", "", rates...], " & ")
        println("$ratio_row \\\\")
        println("$rate_row \\\\")
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println("      \\item \\textit{Notes.} The split ratio is displayed in the first row, with the success rate shown in the second row. Split ratio is calculated as the number of parcels divided by (number of orders × train_test fraction). \$|\\mathcal{I}|\$ represents the number of SKUs and \$|O|\$ represents the number of orders. We used an octa-core AMD 5800X3D CPU with 64 GB RAM.")
println("      \\item \$^a\$ No solution could be found within 3600 seconds.")
println("\\end{tablenotes}")
println("\\end{threeparttable}")

# Print raw data for verification
println("\nRaw split ratio data:")
println(split_ratio_by_combo)
