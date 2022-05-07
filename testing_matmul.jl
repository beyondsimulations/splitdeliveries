
using Distributions
using Random
using DataFrames
using SparseArrays
using StatsPlots

orders         = 50000
min_dependence = 0.01
max_dependence = 0.80
max_group_size = 50
group_link     = 0.05
ind_chance     = 0.20
one_direction  = 0.30
multi_relatio  = 0.50
skus           = 200

transactions = RANDOMTRANS(skus::Int64,
                     orders::Int64,
                     max_group_size::Int64,
                     min_dependence::Float64,
                     max_dependence::Float64,
                     group_link::Float64,
                     ind_chance::Float64,
                     one_direction::Float64,
                     multi_relatio::Float64)
    
sum(transactions)