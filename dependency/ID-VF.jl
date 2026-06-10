# Independent, Variable Frequency (ID-VF)
# No co-purchase dependencies, 100% independent SKUs
# SKU frequency: variable (some SKUs much more popular than others)
include(joinpath(@__DIR__, "defaults.jl"))

sku_frequency_mode = :zipf     # Variable SKU selection frequency
ratio_strong = 0.00      # No strong-dependency groups
ratio_medium = 0.00      # No medium-dependency groups
group_link = 0.00      # No cross-group bridge links
dep_activation_prob = 0.00      # No dependencies to activate
