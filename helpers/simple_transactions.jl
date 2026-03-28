## import packages
include("load_packages.jl")

include("dependency/HD-VF.jl")

X, C, group_sizes = RANDOMTRANS(
    100::Int64,
    1000::Int64,
    mean_order_size::Float64,
    min_order_size::Int64,
    nbd_dispersion::Float64,
    sku_frequency_mode::Symbol,
    zipf_exponent::Float64,
    ceil(Int64, max(100 / group_size_scaling, group_size_min))::Int64,
    ceil(Int64, max(ceil(Int64, max(100 / group_size_scaling, group_size_min)) / mean_group_divisor, mean_group_min))::Int64,
    ratio_strong::Float64,
    ratio_medium::Float64,
    dep_strength_strong::Tuple{Float64,Float64},
    dep_strength_medium::Tuple{Float64,Float64},
    group_link::Float64,
    one_direction::Float64,
    multi_relatio::Float64,
    dep_activation_prob::Float64)


# Convert sparse matrix to DataFrame
df = DataFrame(Matrix(X), :auto)
CSV.write("transactions_sample.csv", df)
