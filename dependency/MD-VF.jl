# Medium Dependencies, Variable Frequency (MD-VF)
# 20% strong groups, 60% medium groups, 20% independent SKUs
# SKU frequency: variable (some SKUs much more popular than others)
include(joinpath(@__DIR__, "defaults.jl"))

sku_frequency_mode = :zipf     # Variable SKU selection frequency
ratio_strong = 0.20      # 20% of SKUs in strong-dependency groups
ratio_medium = 0.60      # 60% of SKUs in medium-dependency groups
group_link = 0.04      # Fraction of cross-group bridge links
dep_activation_prob = 0.60      # Probability that dependencies activate per seed SKU
