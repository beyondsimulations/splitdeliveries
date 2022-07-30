## import packages
    include("load_packages.jl")

# Choose the benchmark which should be evaluated
## Benchmarks used in our article:
### b1_40skus
### b2_100skus
### b3_1000skus
### b4_10000skus
## Dependencies used in our article:
### ID-EF
### MD-EF
### HD-EF
### ID-VF
### MD-VF
### HD-VF
    experiment = "s1_1000skus"
    dependency = "MD-VF"

#  Specify the number of orders and the ratio between test
## and training data for the generated transactional data sets
    train_test = 0.50
    #order_sets  = [round(Int, 1000 * 1/train_test * x) for x =10:10:50]
    order_sets  = [round(Int, 1000 * 1/train_test * x) for x =10:10:100]

# load the data that specifies the dependencies
    include("dependency/$dependency.jl")

# Set the number of cpu cores your computer has at its disposal
    cpu_cores  = 8
    ren_lock = ReentrantLock()

# Choose Optimisations and Heuristics to evaluate in the benchmark
    start = DataFrame(QMKO  = [0], # quadratic-multiple knapsack heuristic with CPLEX as solver
                      QMK   = [1], # quadratic-multiple knapsack heuristic with SBB as solver
                      CHIM  = [1], # main chi-square heuristic without local search
                      CHI   = [1], # chi-square heuristic + local search based on the QMK objective function
                      KL    = [0], # K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2019.07.004
                      KLQ   = [1], # K-LINK optimisation with SBB by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2019.07.004
                      GO    = [1], # greedy orders heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      GP    = [1], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      GS    = [1], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      BS    = [1], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
                      OPT   = [0]) # optimisation model to determine the optimal solution with CPLEX

# Parameters for the KLINK heuristic
## trials: number of different trials with a completly new random solution
## stagnant: number of local search iterations without improvement till termination
## strategy: select the local search strategy of the KLINK heuristic
###          -> 1: move SKU to a different warehouse if it improves the objective and the warehouse has capacity cap_left
###          -> 2: pair-wise exchange between SKUs if it improves the objective
###          -> 3: both strategies above are applied with a chance of 50:50
## klinkstatus: display the current iteration during the heuristic
    trials   = 100
    stagnant = 10
    strategy = 3
    klinkstatus = 0

# Parameters for all Optimisations
## abort: number of seconds until the Optimisation is aborted
## show_opt: specify whether the status of the optimisation should be shown
## allowed_gap: specify the termination criterion in case a gap is allowed in the optimisation
## max_nodes: maximum number of nodes till termination
    abort       = 3600
    show_opt    = false
    allowed_gap = 0.00000
    max_nodes   = 10000000

# Parameters for CHISQUARE
## sig_levels:  significance levels alpha to apply with the chi-square tests
## max_ls:      maximum number of local search runs before termination
## chi_status:  choose whether a detailled progress of the chi heuristic should be shown
    #sig_levels = [1.0e-8]
    sig_levels = [1.0/(10^x) for x = 0:1:9]
    max_ls = 100
    chistatus = false

# Parameters for RANDOM
## iterations: number of different random allocations for the comparison
    iterations = 100

# Parameters for the whole Benchmark
    benchiterations = 20

# Initialise the basic problem by loading the respective capacity constellations
## capacity_benchmark: capacity matrix with column = capacity and row = constellation
    capacity_benchmark  = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark      = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))
    
# Run the benchmark
    print("\n\n### Benchmark of dependency ",dependency," on experiment ",experiment," ###")
    print("\n     benchmark started at ",now(),".\n")
    benchmark = BENCHMARK(capacity_benchmark::Array{Int64,2},
                            skus_benchmark::Vector{Int64},
                            start::DataFrame,
                            order_sets::Vector{Int64},
                            max_dependence::Float64,
                            trials::Int64,
                            stagnant::Int64,
                            strategy::Int64,
                            klinkstatus::Int64,
                            abort::Int64,
                            iterations::Int64,
                            show_opt::Bool,
                            cpu_cores::Int64,
                            allowed_gap::Float64,
                            max_nodes::Int64,
                            sig_levels::Vector{Float64},
                            max_ls::Int64,
                            chistatus::Bool,
                            benchiterations::Int64)
                                                                        
    print("\nbenchmark finished at ",now(),".")
    