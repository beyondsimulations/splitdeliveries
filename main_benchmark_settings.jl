## import packages
include("load_packages.jl")

ren_lock = ReentrantLock()

# Default weight mode (override by setting weight_mode before include)
if !@isdefined(weight_mode)
    weight_mode = :uniform
end

# Defaults for a targeted re-run (override by setting them before include)
## gate_ratio and apply_ls: see their sections below
## results_tag: appended to the names of the result csv and log files
## only_modes: if defined, run only these columns of start, e.g. [:KL]
if !@isdefined(gate_ratio)
    gate_ratio = 2.0
end
if !@isdefined(apply_ls)
    apply_ls = true
end
if !@isdefined(results_tag)
    results_tag = ""
end

# iterate over all dependencies
for dependency in dependencies

    # load the data that specifies the dependencies
    include("dependency/$dependency.jl")

    #  Specify the number of orders and the ratio between test
    ## and training data for the generated transactional data sets
    train_test = 0.50
    order_sets = [2, 10, 20]

    # Set the number of cpu cores your computer has at its disposal
    cpu_cores = 4

    # Choose Optimisations and Heuristics to evaluate in the benchmark
    # Solver-based methods only run at the small scales and KL up to 10,000 SKUs,
    # as it needs a dense SKU-by-SKU link matrix; the scalable heuristics run
    # everywhere. The experiment variable is set by the run scripts generated
    # through run_benchmarks.sh.
    small_scales = experiment in ("100", "1000")
    start = DataFrame(;
        QMK = [0], # quadratic-multiple knapsack heuristic with Gurobi as solver
        QMKJ = [small_scales ? 1 : 0], # quadratic-multiple knapsack heuristic with Juniper as solver
        QMKS = [0], # quadratic-multiple knapsack heuristic with SCIP as solver
        CHI = [1], # main chi-square heuristic
        KL = [experiment in ("100", "1000", "10000") ? 1 : 0], # K-LINK heuristic by S. Zhu, X. Hu, K. Huang and Y. Yuan (2021) https://doi.org/10.1016/j.ejor.2020.08.024
        KLQ = [0], # K-LINK optimisation with Gurobi by S. Zhu, X. Hu, K. Huang and Y. Yuan (2021) https://doi.org/10.1016/j.ejor.2020.08.024
        GO = [1], # greedy orders heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        GP = [1], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        GS = [1], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        BS = [1], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
        EMCI = [1], # extended MCI for D warehouses by Lin et al. (2025)
        IIH = [small_scales ? 1 : 0], # iterative improvement heuristic with Gurobi by Lin et al. (2025) -- 2-warehouse overlapping
        IIHS = [0], # IIH with SCIP by Lin et al. (2025) -- 2-warehouse overlapping
        OPT = [small_scales ? 1 : 0], # optimisation model to determine the optimal solution with Gurobi
    )
    @isdefined(only_modes) && (start[1, Not(only_modes)] .= 0)

    # SKU weight mode (set before include to override)
    ## :uniform  — all weights = 1 (default, classic benchmark)
    ## :frequency — weight = number of orders containing that SKU in training data
    ## :random   — weight sampled uniformly from weight_range
    weight_range = 1:10

    # Parameters for the local search heuristic
    ## apply_ls: apply the local search heuristic after the greedy heuristics (default true, set at the top)

    # Parameters for the KLINK heuristic
    ## trials: maximum number of trials, each starting from a new random solution (Trial_max)
    ## stagnant: stop after this many consecutive trials without a better solution than the best so far (Zhu et al. (2021) name 10 as an example; with stagnant >= trials every trial runs)
    ## strategy: local search of each trial (Fig. 6 of Zhu et al. (2021)), which draws a random warehouse, applies the best improving move of its SKUs and stops when no warehouse has one
    ###          -> 1: move SKU to a different warehouse if it improves the objective, the target has free capacity and the source keeps at least one SKU
    ###          -> 2: pair-wise exchange between SKUs of two warehouses if it improves the objective
    ###          -> 3: the best move of strategy 1 if one improves, otherwise the best move of strategy 2
    ## klinkstatus: display the current iteration during the heuristic
    trials = 100
    stagnant = 100
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
    ## gate_ratio:  factor kappa of the aggregate independence gate (0.0 disables the gate; default 2.0, set at the top)
    sig_levels = [1.0e-5]
    max_ls = 100
    chistatus = false
    min_effect = 0.01
    ls_neighborhood = 1.0

    # Parameters for IIH (Lin et al. 2025)
    ## max_iih_iterations: maximum number of alternating optimization rounds
    ## epsilon_iih: minimum improvement in split deliveries to continue
    max_iih_iterations = 50
    epsilon_iih = 1e-6

    # Parameters for RANDOM
    ## iterations: number of different random allocations for the comparison
    iterations = 10

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
        if !isfile("capacity/capacity_$experiment.csv") ||
            !isfile("capacity/skus_$experiment.csv")
            throw(
                ArgumentError("Capacity or SKUs file not found for experiment: $experiment")
            )
        end
    end

    # Initialise the basic problem by loading the respective capacity constellations
    ## capacity_benchmark: capacity matrix with column = capacity and row = constellation
    capacity_benchmark = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))
    diff_benchmark = vec(readdlm("capacity/diff_$experiment.csv", ';', Float64))
    buff_benchmark = vec(readdlm("capacity/buff_$experiment.csv", ';', Float64))

    # Run the benchmark
    print(
        "\n\n### Benchmark of dependency ",
        dependency,
        " on experiment ",
        experiment,
        " ###",
    )
    print("\n     benchmark started at ", now(), ".\n")
    benchmark = BENCHMARK(
        capacity_benchmark::Array{Int64,2},
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
        min_effect::Float64,
        ls_neighborhood::Float64,
        apply_ls::Bool,
        weight_mode::Symbol,
        weight_range::UnitRange{Int64},
        max_iih_iterations::Int64,
        epsilon_iih::Float64,
        benchiterations::Int64,
        train_test::Float64,
        dependency::String;
        gate_ratio = gate_ratio,
        results_tag = results_tag,
    )

    print("\nbenchmark finished at ", now(), ".")
end
