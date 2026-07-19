using Pkg
Pkg.activate("splitdeliveries")
using CSV, DataFrames, Statistics

# Ablation of the independence gate: CHI with the gate (production
# configuration) against CHI with the gate disabled, on identical scenario
# grids. Uniform storage requirements only.
df_gated = CSV.read("results/overall_results.csv", DataFrame)
df_ungated = CSV.read("results/overall_results_chiungated.csv", DataFrame)

function chi_frame(df)
    out = df[(df.mode .== "CHI_1.0e-5") .& (df.weight_mode .== "uniform"), :]
    out.split_ratio = out.parcel_test ./ out.orders
    out.train_split_ratio = out.parcel_train ./ out.orders
    return out
end

gated = chi_frame(df_gated)
ungated = chi_frame(df_ungated)

# Compare on identical scenario grids only: the re-run of the 14 crashed
# tail scenarios had CHI disabled, so they exist gated-only and would skew
# the per-cell means.
scen = [:dependency, :skus, :wareh, :diff, :buffer, :orders]
common = innerjoin(unique(gated[:, scen]), unique(ungated[:, scen]); on = scen)
gated = innerjoin(gated, common; on = scen)
ungated = innerjoin(ungated, common; on = scen)
println("scenarios compared: $(nrow(common)) (gate-only scenarios dropped)")

# Get unique dependency structures in desired order (ID, MD, HD)
all_deps = unique(gated.dependency)
dependency_levels = []
for prefix in ["ID", "MD", "HD"]
    for dep in sort(all_deps)
        if startswith(dep, prefix)
            push!(dependency_levels, dep)
        end
    end
end

sku_levels = sort(unique(gated.skus))

mean_ratio(frame, dep, sku, col) = begin
    subset = frame[(frame.dependency .== dep) .& (frame.skus .== sku), :]
    nrow(subset) > 0 ? round(mean(subset[!, col]) * 100; digits = 2) : nothing
end

println("Gate ablation analysis (test split ratio, gated vs ungated):")
println("="^60)
for dep in dependency_levels
    for sku in sku_levels
        g = mean_ratio(gated, dep, sku, :split_ratio)
        u = mean_ratio(ungated, dep, sku, :split_ratio)
        gt = mean_ratio(gated, dep, sku, :train_split_ratio)
        ut = mean_ratio(ungated, dep, sku, :train_split_ratio)
        println("$dep $sku: gated $g (train $gt) | ungated $u (train $ut)")
    end
end

# Display the table
println("\n\n")
println(
    "\\caption{Average split ratio (in \\%) of CHI with and without the independence gate}",
)
println("\\label{tab:gate_ablation}")
println("\\begin{threeparttable}")
println("\\begin{tabular}{l" * "rr"^length(sku_levels) * "}")
println("\\toprule")
print(" ")
for sku in sku_levels
    sku_formatted = replace(string(sku), r"(?<=\d)(?=(\d{3})+(?!\d))" => ",")
    print(" & \\multicolumn{2}{c}{$sku_formatted}")
end
println(" \\\\")
print(" ")
for sku in sku_levels
    print(" & gated & ungated")
end
println(" \\\\")
println("\\midrule")

for (i, dep) in enumerate(dependency_levels)
    print("$dep")
    for sku in sku_levels
        g = mean_ratio(gated, dep, sku, :split_ratio)
        u = mean_ratio(ungated, dep, sku, :split_ratio)
        print(" & ", g === nothing ? "-" : "$g")
        print(" & ", u === nothing ? "-" : "$u")
    end
    println(" \\\\")
    if i == 2 || i == 4
        println("\\midrule")
    end
end

println("\\bottomrule")
println("\\end{tabular}")
println("\\begin{tablenotes}")
println("      \\smaller")
println(
    "      \\item \\textit{Notes.} Test split ratios as percentages (lower is better), averaged over all uniform-weight scenarios per dataset structure and SKU count. Column pairs compare CHI with the independence gate (production configuration) against CHI with the gate disabled.",
)
println("\\end{tablenotes}")
println("\\end{threeparttable}")
