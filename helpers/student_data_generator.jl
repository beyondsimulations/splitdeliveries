# Generates the small student practice datasets (20 SKUs, 1000 orders) for
# the three dependency variants under the chosen frequency mode (VF or SF).
# Writes transactions_<V>.csv and dependencies_<V>.csv into student_data/.

const PROJECT_ROOT = joinpath(@__DIR__, "..")
cd(PROJECT_ROOT)
include(joinpath(PROJECT_ROOT, "load_packages.jl"))

using Random

const STUDENT_DIR = joinpath(PROJECT_ROOT, "student_data")
mkpath(STUDENT_DIR)

# Variants to generate. Pass e.g. ARGS = ["SF"] to only produce SF, or leave
# empty to produce both frequency modes.
freq_suffixes = isempty(ARGS) ? ["VF", "SF"] : ARGS
dep_levels    = ["HD", "MD", "ID"]

const STUDENT_SKUS   = 20
const STUDENT_ORDERS = 1000
const STUDENT_SEED   = 20240506   # deterministic so reruns reproduce the data

sku_cols(n) = ["SKU_" * lpad(string(i), 2, '0') for i in 1:n]

for level in dep_levels, suffix in freq_suffixes
    variant = "$(level)-$(suffix)"
    println("Generating $variant ...")

    Random.seed!(STUDENT_SEED)
    include(joinpath(PROJECT_ROOT, "dependency/$(variant).jl"))

    max_gs  = ceil(Int64, max(STUDENT_SKUS / group_size_scaling, group_size_min))
    mean_gs = ceil(Int64, max(max_gs / mean_group_divisor, mean_group_min))

    trans, C, _ = RANDOMTRANS(
        STUDENT_SKUS, STUDENT_ORDERS,
        mean_order_size, min_order_size, nbd_dispersion,
        sku_frequency_mode, zipf_exponent,
        max_gs, mean_gs,
        ratio_strong, ratio_medium,
        dep_strength_strong, dep_strength_medium,
        group_link, one_direction, multi_relatio,
        dep_activation_prob)

    cols = sku_cols(STUDENT_SKUS)

    tx_df = DataFrame(Int.(Matrix(trans)), cols)
    insertcols!(tx_df, 1, :order_id => 1:STUDENT_ORDERS)
    CSV.write(joinpath(STUDENT_DIR, "transactions_$(variant).csv"), tx_df)

    dep_df = DataFrame(Matrix(C), cols)
    insertcols!(dep_df, 1, :added_sku => cols)
    CSV.write(joinpath(STUDENT_DIR, "dependencies_$(variant).csv"), dep_df)
end

println("Done. Files written to $STUDENT_DIR")
