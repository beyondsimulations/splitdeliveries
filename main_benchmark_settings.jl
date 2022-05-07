## import packages
    include("load_packages.jl")

# Choose the benchmark which should be evaluated
## Benchmarks used in our article:
### 1_solution_gap_10
### 2_sensitivity_ind_100
### 2_sensitivity_md_100
### 2_sensitivity_hd_100
### 3_calculation_time_ind
### 3_calculation_time_md
### 3_calculation_time_hd
    experiment = "3_calculation_time_md"

## Alternatively one could specify new transactional data sets and capacity constellations.
## To see how the data has to be specified take a look at the "capacity_***" and "transactions_***" data.

# Set the number of cpu cores your computer has at its disposal
    cpu_cores  = 8
    ren_lock = ReentrantLock()

# Choose Optimisations and Heuristics to evaluate in the benchmark
    start = DataFrame(QMKOPT = [0], # quadratic-multiple knapsack heuristic with CPLEX as solver
                      QMK    = [1], # quadratic-multiple knapsack heuristic with SBB as solver
                      QMKLOC = [0], # quadratic-multiple knapsack heuristic with SBB as solver + local search based on the QMK objective function
                      CHI    = [1], # chi-square heuristic 
                      CHILOC = [1], # chi-square heuristic + local search based on the QMK objective function
                      KLINK  = [0], # K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2019.07.004
                      GP     = [1], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      GPLOC  = [0], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687 + local search
                      GS     = [1], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      GSLOC  = [0], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687 + local search
                      BS     = [1], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      BSLOC  = [0], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687  + local search
                      OPT    = [0], # optimisation model to determine the optimal solution with CPLEX
                      RND    = [1]) # random allocation of SKUs (cannot be deactivated)

# Parameter for the transaction generation if no transactional data is specified under the path 
# "transactions/transactions_$experiment". The benchmark will generate random and independent transactions
# while "order" specifies the number of transactions in each dataset. The parameter max_dependence 
# specifies the maximal strength of the dependecies between products. max-groupsize specifies the 
# maximal group size while group link specifies the ratio of outlinks of each group. Finally, ind_chance
# can be used to determine the chance that instead of group-allocations an independent product is ordered.
    orders         = 200000
    min_dependence = 0.00
    max_dependence = 0.30
    group_link     = 0.05
    ind_chance     = 0.30
    one_direction  = 0.30
    multi_relatio  = 0.50

# Parameters for the KLINK heuristic
## trials: number of different trials with a completly new random solution
## stagnant: number of local search iterations without improvement till termination
## strategy: select the local search strategy of the KLINK heuristic
###          -> 1: move SKU to a different warehouse if it improves the objective and the warehouse has capacity cap_left
###          -> 2: pair-wise exchange between SKUs if it improves the objective
## klinkstatus: display the current iteration during the heuristic
    trials   = 100
    stagnant = 10
    strategy = 2
    klinkstatus = 0

# Parameters for all Optimisations
## abort: number of seconds until the Optimisation is aborted
## show_opt: specify whether the status of the optimisation should be shown
## allowed_gap: specify the termination criterion in case a gap is allowed in the optimisation
## max_nodes: maximum number of nodes till termination
    abort       = 43200
    show_opt    = 0
    allowed_gap = 0.00000
    max_nodes   = 10000000

# Parameters for CHISQUARE
## sig: significance level alpha for the chi-square tests
    sig = 0.05

# Parameters for RANDOM
## iterations: number of different random allocations for the comparison
    iterations = 10

# Initialise the basic problem by loading the respective capacity constellations
## capacity_benchmark: capacity matrix with column = capacity and row = constellation
    capacity_benchmark  = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark      = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))

# Run the benchmark
    parcels_benchmark, 
    time_benchmark, 
    cap_used, 
    parcel_reduction, 
    split_reduction, 
    gap_optimisation = BENCHMARK(capacity_benchmark::Array{Int64,2},
                                 skus_benchmark::Vector{Int64},
                                 start::DataFrame,
                                 orders::Int64,
                                 max_dependence::Float64,
                                 trials::Int64,
                                 stagnant::Int64,
                                 strategy::Int64,
                                 klinkstatus::Int64,
                                 abort::Int64,
                                 iterations::Int64,
                                 show_opt::Int64,
                                 cpu_cores::Int64,
                                 allowed_gap::Float64,
                                 max_nodes::Int64,
                                 sig::Float64)
                                                                         

# Export the results
CSV.write("results/parcels_sent_$experiment.csv", parcels_benchmark)
CSV.write("results/duration_$experiment.csv", time_benchmark)
CSV.write("results/capacity_used_$experiment.csv", cap_used)
CSV.write("results/parcel_reduction_$experiment.csv", parcel_reduction)
CSV.write("results/split_reduction_$experiment.csv", split_reduction)
CSV.write("results/optimisation_gap_$experiment.csv", gap_optimisation)