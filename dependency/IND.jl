# Parameter for the transaction generation if no transactional data is specified under the path 
# "transactions/transactions_$experiment". The benchmark will generate random and independent transactions
# while "order" specifies the number of transactions in each dataset. The parameter max_dependence 
# specifies the maximal strength of the dependecies between products. max-groupsize specifies the 
# maximal group size while group link specifies the ratio of outlinks of each group. Finally, ind_chance
# can be used to determine the chance that instead of group-allocations an independent product is ordered.
    min_dependence = 0.00
    max_dependence = 0.00
    group_link     = 0.00
    ind_chance     = 0.00
    one_direction  = 0.00
    multi_relatio  = 0.00