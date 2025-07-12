## import packages
include("load_packages.jl")

include("dependency/HD-VF.jl")

X = RANDOMTRANS(
    100::Int64,
    1000::Int64,
    skus_in_order::Float64,
    sku_frequency::Float64,
    10::Int64,
    min_dependence::Float64,
    max_dependence::Float64,
    group_link::Float64,
    ind_chance::Float64,
    one_direction::Float64,
    multi_relatio::Float64)


# Convert sparse matrix to DataFrame
df = DataFrame(Matrix(X), :auto)
CSV.write("transactions_sample.csv", df)

