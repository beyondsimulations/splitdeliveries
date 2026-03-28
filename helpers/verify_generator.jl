## Verification script for the RANDOMTRANS transaction data generator
## Generates small instances (100 SKUs, 1000 orders) for each of the 6 configs
## and produces diagnostic plots saved to results/verification/

const PROJECT_ROOT = joinpath(@__DIR__, "..")
cd(PROJECT_ROOT)
include(joinpath(PROJECT_ROOT, "load_packages.jl"))

# Create output directory
mkpath("results/verification")

configs = ["HD-SF", "HD-VF", "MD-SF", "MD-VF", "ID-SF", "ID-VF"]

for config in configs
    println("\n", "="^60)
    println("Generating and verifying: $config")
    println("="^60)

    # Load dependency parameters
    include(joinpath(PROJECT_ROOT, "dependency/$config.jl"))

    skus = 100
    orders = 1000

    # Generate transactions
    max_gs = ceil(Int64, max(skus / group_size_scaling, group_size_min))
    mean_gs = ceil(Int64, max(max_gs / mean_group_divisor, mean_group_min))
    trans, C, group_sizes_gen = RANDOMTRANS(
        skus, orders,
        mean_order_size, min_order_size, nbd_dispersion,
        sku_frequency_mode, zipf_exponent,
        max_gs, mean_gs,
        ratio_strong, ratio_medium,
        dep_strength_strong, dep_strength_medium,
        group_link, one_direction, multi_relatio,
        dep_activation_prob)

    # Compute basic statistics
    order_sizes = vec(sum(trans, dims=2))
    sku_freqs = vec(sum(trans, dims=1))

    # --- Summary statistics ---
    println("\nConfig: $config")
    println("  Orders: $(size(trans,1)), SKUs: $(size(trans,2))")
    println("  Mean order size: $(round(mean(order_sizes), digits=2)) (target: $mean_order_size)")
    println("  Min order size: $(minimum(order_sizes)) (target: $min_order_size)")
    println("  Max order size: $(maximum(order_sizes))")
    println("  Median order size: $(median(order_sizes))")
    println("  SKU freq std/mean (CV): $(round(std(sku_freqs)/mean(sku_freqs), digits=3))")
    println("  Non-zero entries in C: $(nnz(C))")
    if length(group_sizes_gen) > 0
        println("  Groups built: $(length(group_sizes_gen))")
        println("  Mean group size: $(round(mean(group_sizes_gen), digits=2))")
        println("  Max group size: $(maximum(group_sizes_gen))")
    else
        println("  Groups built: 0 (independent config)")
    end

    # --- Plot 1: Order size histogram ---
    p1 = histogram(order_sizes,
        title="Order Size ($config)",
        xlabel="Items per order", ylabel="Frequency",
        legend=false, bins=maximum(order_sizes) - min_order_size + 1,
        color=:steelblue, linecolor=:white)
    vline!(p1, [mean_order_size], color=:red, linestyle=:dash, linewidth=2,
        label="Target mean")
    savefig(p1, "results/verification/$(config)_order_sizes.png")

    # --- Plot 2: SKU frequency (sorted descending) ---
    p2 = bar(sort(sku_freqs, rev=true),
        title="SKU Frequency ($config)",
        xlabel="SKU (sorted by frequency)", ylabel="Appearances",
        legend=false, color=:darkorange, linecolor=:white)
    savefig(p2, "results/verification/$(config)_sku_frequency.png")

    # --- Plot 3: Dependency matrix C heatmap ---
    p3 = heatmap(Matrix(C),
        title="Dependency Matrix C ($config)",
        xlabel="SKU", ylabel="SKU", color=:viridis,
        aspect_ratio=:equal, size=(600,500))
    savefig(p3, "results/verification/$(config)_dependency_matrix.png")

    # --- Plot 4: Empirical co-purchase correlation matrix ---
    trans_float = Float64.(Matrix(trans))
    cooccurrence = trans_float' * trans_float
    diag_vals = sqrt.(diag(cooccurrence) .+ 1.0)
    D_inv = Diagonal(1.0 ./ diag_vals)
    corr = D_inv * cooccurrence * D_inv
    p4 = heatmap(corr,
        title="Co-purchase Correlation ($config)",
        xlabel="SKU", ylabel="SKU", color=:viridis,
        aspect_ratio=:equal, size=(600,500))
    savefig(p4, "results/verification/$(config)_copurchase_correlation.png")

    # --- Plot 5: Group size distribution ---
    if length(group_sizes_gen) > 0
        p5 = histogram(group_sizes_gen,
            title="Group Size Distribution ($config)",
            xlabel="Group size", ylabel="Count",
            legend=false, color=:seagreen, linecolor=:white)
        savefig(p5, "results/verification/$(config)_group_sizes.png")
    end

    # --- Plot 6: Degree distribution in C ---
    degrees = vec(sum(C .> 0, dims=1))
    if any(degrees .> 0)
        p6 = histogram(degrees[degrees .> 0],
            title="SKU Out-Degree in C ($config)",
            xlabel="Number of connections", ylabel="Count",
            legend=false, color=:mediumpurple, linecolor=:white)
        savefig(p6, "results/verification/$(config)_degree_distribution.png")
    end

    # --- Combined summary plot ---
    plots_to_combine = [p1, p2, p3, p4]
    p_combined = plot(plots_to_combine..., layout=(2,2), size=(1200,1000),
        plot_title="Transaction Generator Verification: $config")
    savefig(p_combined, "results/verification/$(config)_summary.png")

    println("  Plots saved to results/verification/$(config)_*.png")
end

println("\n", "="^60)
println("Verification complete. Check results/verification/ for plots.")
println("="^60)
