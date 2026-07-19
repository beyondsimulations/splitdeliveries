# Targeted re-run of the 14 uniform 1,000-SKU scenarios missing from run 2
# (the K = 8/10 buffered tail; constellations 17-20 of capacity_1000.csv).
# CHI stays OFF: its rows come from the gated re-run in results_chigate/ and
# are swapped in by aggregate_benchmarks.jl. Everything else mirrors
# main_benchmark_settings.jl for experiment "1000" with uniform weights.
#
# Run from the repo root on the benchmark machine:
#   julia run_missing_tail.jl
#
# The run writes its own files (results/1000fix_benchmark_*.csv, kept as a
# record) and appends the new rows to the canonical
# results/1000_benchmark_*.csv (a .bak backup is created first). Afterwards
# re-run: julia aggregate_benchmarks.jl

include("load_packages.jl")

# (dependency, wareh, diff, buffer, order_set) of every missing scenario,
# as identified from overall_results.csv on 2026-07-17
missing_scenarios = [
    ("HD-VF", 10, 0.0, 0.2, 20),
    ("HD-VF", 10, 0.2, 0.2, 2),
    ("HD-VF", 10, 0.2, 0.2, 10),
    ("HD-VF", 10, 0.2, 0.2, 20),
    ("ID-VF", 8, 0.0, 0.2, 20),
    ("ID-VF", 8, 0.2, 0.2, 2),
    ("ID-VF", 8, 0.2, 0.2, 10),
    ("ID-VF", 8, 0.2, 0.2, 20),
    ("ID-VF", 10, 0.0, 0.2, 2),
    ("ID-VF", 10, 0.0, 0.2, 10),
    ("ID-VF", 10, 0.0, 0.2, 20),
    ("ID-VF", 10, 0.2, 0.2, 2),
    ("ID-VF", 10, 0.2, 0.2, 10),
    ("ID-VF", 10, 0.2, 0.2, 20),
]

# Filenames carry "1000fix"; the instance seeds depend on the ORIGINAL
# constellation row indices, which only_pairs preserves
experiment = "1000fix"
weight_mode = :uniform

capacity_benchmark = readdlm("capacity/capacity_1000.csv", ';', Int64)
skus_benchmark = vec(readdlm("capacity/skus_1000.csv", ';', Int64))
diff_benchmark = vec(readdlm("capacity/diff_1000.csv", ';', Float64))
buff_benchmark = vec(readdlm("capacity/buff_1000.csv", ';', Float64))

wareh_of(a) = count(x -> x > 0, capacity_benchmark[a, :])

for dependency in unique(first.(missing_scenarios))
    include("dependency/$dependency.jl")

    only_pairs = Set{Tuple{Int,Int}}()
    for (dep, K, dif, buf, os) in missing_scenarios
        dep == dependency || continue
        rows = [
            a for a in axes(capacity_benchmark, 1) if
            wareh_of(a) == K && diff_benchmark[a] == dif && buff_benchmark[a] == buf
        ]
        length(rows) == 1 ||
            error("expected one constellation for K=$K diff=$dif buff=$buf, got $rows")
        push!(only_pairs, (rows[1], os))
    end
    print("\n### $dependency: re-running pairs $(sort(collect(only_pairs))) ###\n")

    train_test = 0.50
    order_sets = [2, 10, 20]
    cpu_cores = 4
    start = DataFrame(;
        QMK = [0], QMKJ = [1], QMKS = [0], CHI = [0], KL = [1], KLQ = [0],
        GO = [1], GP = [1], GS = [1], BS = [1], EMCI = [1],
        IIH = [0], IIHS = [0], OPT = [1],
    )
    weight_range = 1:10
    apply_ls = true
    trials = 100
    stagnant = 10
    strategy = 3
    klinkstatus = false
    abort = 900
    show_opt = false
    allowed_gap = 0.00000
    max_nodes = 10000000
    sig_levels = [1.0e-5]
    max_ls = 100
    chistatus = false
    min_effect = 0.01
    ls_neighborhood = 1.0
    gate_ratio = 2.0
    max_iih_iterations = 50
    epsilon_iih = 1e-6
    iterations = 10
    benchiterations = 1

    fixed = BENCHMARK(
        capacity_benchmark,
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
        min_effect,
        ls_neighborhood,
        apply_ls,
        weight_mode,
        weight_range,
        max_iih_iterations,
        epsilon_iih,
        benchiterations,
        train_test,
        dependency;
        gate_ratio = gate_ratio,
        only_pairs = only_pairs,
    )

    # Append the new rows to the canonical results file (backup kept)
    canonical = "results/1000_benchmark_$(dependency).csv"
    cp(canonical, canonical * ".bak"; force = true)
    existing = CSV.read(canonical, DataFrame; stringtype = String)
    CSV.write(canonical, vcat(existing, fixed))
    print("\nAppended $(nrow(fixed)) rows to $canonical (backup: $canonical.bak)\n")
end

print("\n### All missing scenarios re-run. Now run: julia aggregate_benchmarks.jl ###\n")
