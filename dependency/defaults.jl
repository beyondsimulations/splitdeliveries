# Shared default parameters for transaction generation
# Override any of these in the specific config file after including this file

# Order parameters
mean_order_size      = 3.0       # Expected number of unique items per order
min_order_size       = 2         # Minimum items per order (single-item orders excluded)
nbd_dispersion       = 2.0       # NBD dispersion r (r=1 → geometric; r>1 → heavier tail)

# SKU frequency distribution
zipf_exponent        = 0.9       # Zipf exponent α: P(rank k) ∝ 1/k^α

# Dependency strength ranges
dep_strength_strong  = (0.50, 0.90)  # Co-purchase probability range for strong groups
dep_strength_medium  = (0.10, 0.50)  # Co-purchase probability range for medium groups

# Group size scaling
group_size_scaling   = 20        # Max group size divisor: g_max = ⌈max(S / this, group_size_min)⌉
group_size_min       = 10        # Minimum max group size
mean_group_divisor   = 4         # Mean group size divisor: g_mean = ⌈max(g_max / this, mean_group_min)⌉
mean_group_min       = 3         # Minimum mean group size

# Group topology
one_direction        = 0.60      # Probability of asymmetric (one-way) dependency
multi_relatio        = 0.40      # Peripheral-to-peripheral edge density within groups
