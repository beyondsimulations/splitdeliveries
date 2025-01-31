# Transaction Generation Parameters
# These parameters control the generation of random transactions when no transactional data 
# is specified under "transactions/transactions_$experiment"

# Order-Level Parameters
skus_in_order  = 3.00    # Mean of normal distribution for number of SKUs per order
                         # Actual: 1 + floor(abs(rand(Normal(0,skus_in_order))))
sku_frequency  = 3.50    # Controls SKU selection probability distribution
                         # 0.00: Uniform random selection
                         # >0.0: Normal distribution: 1 + floor(Int64, abs(rand(Normal(0,skus/sku_frequency))))

# Dependency Parameters
min_dependence = 0.00    # Minimum correlation strength between related SKUs
max_dependence = 0.20    # Maximum correlation strength between related SKUs

# Group Structure Parameters
group_link     = 0.02    # Ratio of SKUs that get random connections outside their group
ind_chance     = 0.30    # Ratio of SKUs not assigned to any cluster
one_direction  = 0.80    # Probability of one-way (vs two-way) dependencies between SKUs
multi_relatio  = 0.20    # Depth of relations within product groups (controls multi-dimensional dependencies)