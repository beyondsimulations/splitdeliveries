# Quick test benchmark: 100 SKUs, all 6 dependencies, no Gurobi-dependent heuristics
# Run from project root: julia helpers/test_benchmark.jl

const PROJECT_ROOT = joinpath(@__DIR__, "..")
cd(PROJECT_ROOT)

# Set experiment and dependencies before loading the benchmark
experiment = "1000"
dependencies = ["HD-SF", "HD-VF", "MD-SF", "MD-VF", "ID-SF", "ID-VF"]

include(joinpath(PROJECT_ROOT, "load_packages.jl"))

ren_lock = ReentrantLock()

for dependency in dependencies

    include(joinpath(PROJECT_ROOT, "dependency/$dependency.jl"))

    train_test = 0.50
    order_sets = [4]          # only 1 order set for speed
    cpu_cores = 4
    benchiterations = 1

    # Disable all Gurobi-dependent heuristics
    start = DataFrame(
        QMK=[0],    # requires Gurobi
        QMKJ=[1],   # Juniper (no license needed)
        QMKS=[0],   # SCIP
        CHI=[1],    # chi-square heuristic
        KL=[1],     # K-LINK heuristic
        KLQ=[0],    # K-LINK with Gurobi
        GO=[1],     # greedy orders
        GP=[1],     # greedy pairs
        GS=[1],     # greedy seeds
        BS=[1],     # bestselling
        EMCI=[1],   # extended MCI
        IIH=[0],    # requires Gurobi
        IIHS=[0],   # requires SCIP
        OPT=[0],    # exact optimization (Gurobi)
    )

    trials = 10
    stagnant = 5
    strategy = 3
    klinkstatus = false

    abort = 60
    show_opt = false
    allowed_gap = 0.00
    max_nodes = 100000

    sig_levels = [1.0e-2]
    max_ls = 10
    chistatus = false

    max_iih_iterations = 10
    epsilon_iih = 1e-6

    iterations = 10

    capacity_benchmark = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))
    diff_benchmark = vec(readdlm("capacity/diff_$experiment.csv", ';', Float64))
    buff_benchmark = vec(readdlm("capacity/buff_$experiment.csv", ';', Float64))

    print("\n\n### Test Benchmark: $dependency on experiment $experiment ###")
    print("\n     started at ", now(), ".\n")

    benchmark = BENCHMARK(capacity_benchmark,
        skus_benchmark,
        diff_benchmark,
        buff_benchmark,
        start,
        order_sets,
        trials,
        stagnant,
        strategy,
        klinkstatus,
        abort,
        iterations,
        show_opt,
        cpu_cores,
        allowed_gap,
        max_nodes,
        sig_levels,
        max_ls,
        chistatus,
        max_iih_iterations,
        epsilon_iih,
        benchiterations,
        train_test,
        dependency)

    print("\n$dependency finished at ", now(), ".\n")
end

println("\n\nAll test benchmarks complete!")
