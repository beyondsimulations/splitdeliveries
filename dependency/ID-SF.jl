# Transaction Generation Parameters — Independent, Static Frequency (ID-SF)
# Dependency structure: fully independent (no co-purchase dependencies)
# SKU frequency: uniform (all SKUs equally likely as seed)

# Order parameters
mean_order_size      = 2.5       # Expected number of unique items per order
min_order_size       = 2         # Minimum items per order (single-item orders excluded)
nbd_dispersion       = 1.0       # NBD dispersion r (r=1 → geometric; r>1 → heavier tail)

# SKU frequency distribution
sku_frequency_mode   = :uniform  # Uniform random SKU selection
zipf_exponent        = 1.0       # Unused when mode = :uniform

# Dependency structure — fully independent
ratio_strong         = 0.00      # No strong-dependency groups
ratio_medium         = 0.00      # No medium-dependency groups (100% independent)
dep_strength_strong  = (0.0, 0.0)
dep_strength_medium  = (0.0, 0.0)
group_size_scaling   = 20        # Max group size divisor (unused for ID)
group_size_min       = 10        # Minimum max group size (unused for ID)
mean_group_divisor   = 4         # Mean group size divisor (unused for ID)
mean_group_min       = 3         # Minimum mean group size (unused for ID)
group_link           = 0.00      # No cross-group links
one_direction        = 0.00      # Unused
multi_relatio        = 0.00      # Unused
dep_activation_prob  = 0.00      # Unused (no dependencies to activate)
