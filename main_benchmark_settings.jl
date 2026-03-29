## import packages
include("load_packages.jl")

ren_lock = ReentrantLock()

# iterate over all dependencies
for dependency in dependencies

    # load the data that specifies the dependencies
    include("dependency/$dependency.jl")

    #  Specify the number of orders and the ratio between test
    ## and training data for the generated transactional data sets
    train_test = 0.50
    order_sets = [4, 8, 12]

    # Set the number of cpu cores your computer has at its disposal
    cpu_cores = 4

    # Choose Optimisations and Heuristics to evaluate in the benchmark
    start = DataFrame(
        QMK=[0], # quadratic-multiple knapsack heuristic with Gurobi as solver
        QMKJ=[0], # quadratic-multiple knapsack heuristic with Juniper as solver
        QMKS=[0], # quadratic-multiple knapsack heuristic with SCIP as solver
        CHIM=[1], # main chi-square heuristic without local search
        CHI=[0], # chi-square heuristic + local search based on the QMK objective function
        KL=[0], # K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2020.08.024
        KLQ=[0], # K-LINK optimisation with Gurobi by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2020.08.024
        GO=[0], # greedy orders heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        GP=[0], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        GS=[0], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        BS=[0], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        EMCI=[1], # extended MCI for D warehouses by Lin et al. (2025)
        IIH=[0], # iterative improvement heuristic with Gurobi by Lin et al. (2025) -- 2-warehouse overlapping
        IIHS=[0], # IIH with SCIP by Lin et al. (2025) -- 2-warehouse overlapping
        OPT=[0], # optimisation model to determine the optimal solution with Gurobi
    )

    # Parameters for the KLINK heuristic
    ## trials: number of different trials with a completly new random solution
    ## stagnant: number of local search iterations without improvement till termination
    ## strategy: select the local search strategy of the KLINK heuristic
    ###          -> 1: move SKU to a different warehouse if it improves the objective and the warehouse has capacity cap_left
    ###          -> 2: pair-wise exchange between SKUs if it improves the objective
    ###          -> 3: both strategies above are applied with a chance of 50:50
    ## klinkstatus: display the current iteration during the heuristic
    trials = 100
    stagnant = 10
    strategy = 3
    klinkstatus = false

    # Parameters for all Optimisations
    ## abort: number of seconds until the Optimisation is aborted
    ## show_opt: specify whether the status of the optimisation should be shown
    ## allowed_gap: specify the termination criterion in case a gap is allowed in the optimisation
    ## max_nodes: maximum number of nodes till termination
    abort = 900
    show_opt = false
    allowed_gap = 0.00000
    max_nodes = 10000000

    # Parameters for CHISQUARE
    ## sig_levels:  significance levels alpha to apply with the chi-square tests
    ## max_ls:      maximum number of local search runs before termination
    ## chi_status:  choose whether a detailled progress of the chi heuristic should be shown
    sig_levels = [1.0e-2]
    max_ls = 100
    chistatus = false

    # Parameters for IIH (Lin et al. 2025)
    ## max_iih_iterations: maximum number of alternating optimization rounds
    ## epsilon_iih: minimum improvement in split deliveries to continue
    max_iih_iterations = 50
    epsilon_iih = 1e-6

    # Parameters for RANDOM
    ## iterations: number of different random allocations for the comparison
    iterations = 100

    # Parameters for the whole Benchmark
    benchiterations = 1

    # Configuration validation
    function validate_settings()
        if !all(0 .<= sig_levels .<= 1)
            throw(ArgumentError("Significance levels must be between 0 and 1"))
        end
        if cpu_cores < 1
            throw(ArgumentError("CPU cores must be positive"))
        end
        if !isfile("capacity/capacity_$experiment.csv") || !isfile("capacity/skus_$experiment.csv")
            throw(ArgumentError("Capacity or SKUs file not found for experiment: $experiment"))
        end
    end

    # Initialise the basic problem by loading the respective capacity constellations
    ## capacity_benchmark: capacity matrix with column = capacity and row = constellation
    capacity_benchmark = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))
    diff_benchmark = vec(readdlm("capacity/diff_$experiment.csv", ';', Float64))
    buff_benchmark = vec(readdlm("capacity/buff_$experiment.csv", ';', Float64))

    # Run the benchmark
    print("\n\n### Benchmark of dependency ", dependency, " on experiment ", experiment, " ###")
    print("\n     benchmark started at ", now(), ".\n")
    benchmark = BENCHMARK(capacity_benchmark::Array{Int64,2},
        skus_benchmark::Vector{Int64},
        diff_benchmark::Vector{Float64},
        buff_benchmark::Vector{Float64},
        start::DataFrame,
        order_sets::Vector{Int64},
        trials::Int64,
        stagnant::Int64,
        strategy::Int64,
        klinkstatus::Bool,
        abort::Int64,
        iterations::Int64,
        show_opt::Bool,
        cpu_cores::Int64,
        allowed_gap::Float64,
        max_nodes::Int64,
        sig_levels::Vector{Float64},
        max_ls::Int64,
        chistatus::Bool,
        max_iih_iterations::Int64,
        epsilon_iih::Float64,
        benchiterations::Int64,
        train_test::Float64,
        dependency::String)

    print("\nbenchmark finished at ", now(), ".")
end
