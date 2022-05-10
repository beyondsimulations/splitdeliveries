## import packages
    include("load_packages.jl")

# Choose the benchmark which should be evaluated
## Benchmarks used in our article:
### 1_gap
### 2_sensitivity
### 3_time
## Dependencies used in our article:
### IND
### MD
### HD
    experiment = "2_sensitivity"
    dependency = "MD"

#  Specify the number of orders and the ratio between test
## and training data for the generated transactional data sets
    orders     = 500000
    train_test = 0.80

# load the data that specifies the dependencies
    include("dependency/$dependency.jl")

# Set the number of cpu cores your computer has at its disposal
    cpu_cores  = 8
    ren_lock = ReentrantLock()

# Choose Optimisations and Heuristics to evaluate in the benchmark
    start = DataFrame(QMKOPT = [0], # quadratic-multiple knapsack heuristic with CPLEX as solver
                      QMK    = [1], # quadratic-multiple knapsack heuristic with SBB as solver
                      CHI    = [1], # chi-square heuristic 
                      CHILOC = [1], # chi-square heuristic + local search based on the QMK objective function
                      KLINK  = [0], # K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2019.07.004
                      GP     = [1], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      GS     = [1], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      BS     = [1], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      OPT    = [0], # optimisation model to determine the optimal solution with CPLEX
                      RND    = [1]) # random allocation of SKUs (cannot be deactivated)

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
    abort       = 3600
    show_opt    = 0
    allowed_gap = 0.00000
    max_nodes   = 10000000

# Parameters for CHISQUARE
## sig: significance level alpha for the chi-square tests
    sig = 0.01

# Parameters for RANDOM
## iterations: number of different random allocations for the comparison
    iterations = 10

# Initialise the basic problem by loading the respective capacity constellations
## capacity_benchmark: capacity matrix with column = capacity and row = constellation
    capacity_benchmark  = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark      = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))
    
# Run the benchmark
    print("\n\n### Benchmark of dependency ",dependency," on experiment ",experiment," ###")
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
    CSV.write("results/$(experiment)_a_parcels_sent_$dependency.csv",       parcels_benchmark)
    CSV.write("results/$(experiment)_b_duration_$dependency.csv",           time_benchmark)
    CSV.write("results/$(experiment)_c_capacity_used_$dependency.csv",      cap_used)
    CSV.write("results/$(experiment)_d_parcel_reduction_$dependency.csv",   parcel_reduction)
    CSV.write("results/$(experiment)_e_split_reduction_$dependency.csv",    split_reduction)
    CSV.write("results/$(experiment)_f_optimisation_gap_$dependency.csv",   gap_optimisation)
    