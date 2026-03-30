# Medium Dependencies, Static Frequency (MD-SF)
# 20% strong groups, 50% medium groups, 30% independent SKUs
# SKU frequency: uniform (all SKUs equally likely as seed)
include(joinpath(@__DIR__, "defaults.jl"))

sku_frequency_mode   = :uniform  # Uniform random SKU selection
ratio_strong         = 0.20      # 20% of SKUs in strong-dependency groups
ratio_medium         = 0.60      # 50% of SKUs in medium-dependency groups
group_link           = 0.04      # Fraction of cross-group bridge links
dep_activation_prob  = 0.60      # Probability that dependencies activate per seed SKU
