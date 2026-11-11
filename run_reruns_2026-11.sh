#!/bin/bash
# Usage: ./run_reruns_2026-11.sh kl|chi (from the repository root, uniform weights only, needs Julia 1.12.5 via juliaup)
# kl:  KL alone on 100, 1000 and 10000 SKUs -> results/<experiment>_benchmark_<dependency>_klfaithful.csv
# chi: CHI alone without gate and local search on all scales -> results/<experiment>_benchmark_<dependency>_gateoff.csv

GROUP="$1"
if [ "$GROUP" != "kl" ] && [ "$GROUP" != "chi" ]; then
    echo "Usage: ./run_reruns_2026-11.sh kl|chi"
    exit 1
fi

# Folder for the generated per-dependency run scripts
mkdir -p runs
PROJDIR="$(pwd)"

# Array of dependencies
dependencies=(
    "ID-SF"
    "ID-VF"
    "MD-SF"
    "MD-VF"
    "HD-SF"
    "HD-VF"
)

# launch <name> <experiment> <settings>: one tmux session with one window per dependency
launch() {
    session_name="rerun_$1_$2"
    tmux new-session -d -s $session_name

    for dep in "${dependencies[@]}"; do
        # Create a new window
        tmux new-window -t $session_name -n "$dep"

        # Generate run script in runs/ folder
        cat > "runs/rerun_$1_$2_${dep}.jl" << EOL
    include("${PROJDIR}/load_packages.jl")
    dependencies=["${dep}"]
    experiment = "$2"
    weight_mode = :uniform
$3
    include("${PROJDIR}/main_benchmark_settings.jl")
EOL

        tmux send-keys -t $session_name:$dep "cd ${PROJDIR} && julia +1.12.5 runs/rerun_$1_$2_${dep}.jl" C-m
    done

    echo "Started $1 experiment=$2 in tmux session '${session_name}'"
}

if [ "$GROUP" = "kl" ]; then
    for experiment in "100" "1000" "10000"; do
        launch kl $experiment '    only_modes = [:KL]
    results_tag = "_klfaithful"'
    done
else
    for experiment in "100" "1000" "10000" "100000"; do
        launch chi $experiment '    only_modes = [:CHI]
    gate_ratio = 0.0
    apply_ls = false
    results_tag = "_gateoff"'
    done
fi

echo "All re-runs launched. Attach with: tmux attach -t rerun_${GROUP}_100"
