# Transaction Generation Parameters
# These parameters control the generation of random transactions when no transactional data
# is specified under "transactions/transactions_$experiment"

# Order-Level Parameters
skus_in_order  = 3.00   # This parameter sets the standard deviation of the
                        # normal distribution used to generate the number of unique items (SKUs)
                        # in each order. The actual number of SKUs per order is calculated as:
                        # 1 + floor(abs(random_normal(mean=0, std=skus_in_order)))
                        # This ensures orders have at least 1 SKU, with a distribution
                        # that increases in spread as this parameter increases.
sku_frequency  = 3.00    # Controls SKU selection probability distribution
                         # 0.00: Uniform random selection
                         # >0.0: Normal distribution: 1 + floor(Int64, abs(rand(Normal(0,skus/sku_frequency))))

# Dependency Parameters
min_dependence = 0.00    # Minimum correlation strength between related SKUs
max_dependence = 0.00    # Maximum correlation strength between related SKUs

# Group Structure Parameters
group_link     = 0.00    # Ratio of SKUs that get random connections outside their group
ind_chance     = 1.00    # Ratio of SKUs not assigned to any cluster
one_direction  = 0.00    # Probability of one-way (vs two-way) dependencies between SKUs
multi_relatio  = 0.00    # Depth of relations within product groups (controls multi-dimensional dependencies)
