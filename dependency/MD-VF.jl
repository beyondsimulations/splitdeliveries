# Transaction Generation Parameters — Medium Dependencies, Variable Frequency (MD-VF)
# Dependency structure: 20% strong groups, 40% medium groups, 40% independent SKUs
# SKU frequency: Zipf/power-law (some SKUs much more popular than others)

# Order parameters
mean_order_size      = 2.5       # Expected number of unique items per order
min_order_size       = 2         # Minimum items per order (single-item orders excluded)
nbd_dispersion       = 1.0       # NBD dispersion r (r=1 → geometric; r>1 → heavier tail)

# SKU frequency distribution
sku_frequency_mode   = :zipf     # Zipf/power-law SKU selection
zipf_exponent        = 1.0       # Zipf exponent α: P(rank k) ∝ 1/k^α

# Dependency structure (degree-corrected SBM)
ratio_strong         = 0.20      # 20% of SKUs in strong-dependency groups
ratio_medium         = 0.40      # 40% of SKUs in medium-dependency groups (40% independent)
dep_strength_strong  = (0.50, 0.80)  # Co-purchase probability range for strong groups
dep_strength_medium  = (0.20, 0.50)  # Co-purchase probability range for medium groups
group_size_scaling   = 20        # Max group size divisor: g_max = ⌈max(S / this, group_size_min)⌉
group_size_min       = 10        # Minimum max group size
mean_group_divisor   = 4         # Mean group size divisor: g_mean = ⌈max(g_max / this, mean_group_min)⌉
mean_group_min       = 3         # Minimum mean group size
group_link           = 0.03      # Fraction of cross-group bridge links
one_direction        = 0.60      # Probability of asymmetric (one-way) dependency
multi_relatio        = 0.40      # Peripheral-to-peripheral edge density within groups
dep_activation_prob  = 0.60      # Probability that dependencies activate per seed SKU
