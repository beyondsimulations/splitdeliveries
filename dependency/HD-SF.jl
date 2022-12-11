# Parameter for the transaction generation if no transactional data is specified under the path 
# "transactions/transactions_$experiment". The benchmark will generate random transactions. The
# generation is subject to the following parameters for the function "RANDOMTRANS":
## skus:            number of SKUs (products) of the transactional data sets
## orders:          number of orders in the transactional data sets
## max_group_size:  maximal size of a cluster of SKUs with correlations - the actual group sizes
##                  are then drawn from an uniform distribution between 1 and max_group_size
## min_dependence:  specifies the minimal dependence between correlated SKUs of a group
## max_dependence:  specifies the maximal dependence between correlated SKUs of a group - the actual
##                  dependency is then drawn from an uniform distribution between min_dependence and
##                  max_dependence
## group_link:      specifies the ratio of SKUs to all SKUs to which randomly another SKU from the sets
##                  of all SKUs is assigned to prevent "closed" groups. The strength of the dependency 
##                  is also drawn from an uniform distribution between min_dependence and max_dependence
## ind_chance:      the ratio of SKUs not assigned to a cluster of SKUs via max_group_size
## one_direction:   for each pair of SKUs with dependencys the ratio of one_direction decides whether those
##                  dependencies are specified in one direction or in both directions. If a random number between
##                  0 and 1 is higher than one_direction, the relation is true for both directions
## multi_relatio:   for the SKUs in a product group the parameter multi_relatio specifies the depth of the relations
##                  between the SKUs dependening on the overall number of SKUs in the corresponding group_size. Thus,
##                  the dependencies are not only pair-wise but multi-dimensional
## skus_in_order:   the number of SKUs per order is currently drawn from a normal distribution with the following
##                  implementation for each order: skus_order = 1 + floor(abs(rand(Normal(0,skus_in_order))))
## sku_frequency:   the frequency of an SKU appearing in an order independently is currently determined by
##                  the following implementation: new_sku = 1 + floor(Int64, abs(rand(Normal(0,skus/sku_frequency)))).
##                  The exception is sku_frequency = 0 -> in this case new_sku = rand(1:skus).
    skus_in_order  = 3.00
    sku_frequency  = 0.00
    min_dependence = 0.00
    max_dependence = 0.50
    group_link     = 0.04
    ind_chance     = 0.15
    one_direction  = 0.60
    multi_relatio  = 0.40