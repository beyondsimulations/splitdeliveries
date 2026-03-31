#!/bin/bash

# ./run_benchmarks.sh [weight_mode]
# Runs benchmarks in parallel via tmux.
# Usage:
#   ./run_benchmarks.sh              # uniform weights (default)
#   ./run_benchmarks.sh frequency    # frequency-based weights
#   ./run_benchmarks.sh random       # random weights

# Weight mode (default: uniform)
WEIGHT_MODE="${1:-uniform}"

# Array of dependencies
dependencies=(
    "ID-SF"
    "ID-VF"
    "MD-SF"
    "MD-VF"
    "HD-SF"
    "HD-VF"
)

# Experiments to run in parallel
experiments=("100")

for experiment in "${experiments[@]}"; do
    session_name="benchmarks_${experiment}_${WEIGHT_MODE}"
    tmux new-session -d -s $session_name

    for i in "${!dependencies[@]}"; do
        dep="${dependencies[$i]}"

        # Create a new window
        tmux new-window -t $session_name -n "$dep"

        # Generate run script in runs/ folder
        PROJDIR="$(pwd)"
        cat > "runs/run_${experiment}_${dep}_${WEIGHT_MODE}.jl" << EOL
    include("${PROJDIR}/load_packages.jl")
    dependencies=["${dep}"]
    experiment = "${experiment}"
    weight_mode = :${WEIGHT_MODE}
    include("${PROJDIR}/main_benchmark_settings.jl")
EOL

        tmux send-keys -t $session_name:$dep "cd $(pwd) && julia runs/run_${experiment}_${dep}_${WEIGHT_MODE}.jl" C-m
    done

    echo "Started experiment=${experiment} weight_mode=${WEIGHT_MODE} in tmux session '${session_name}'"
done

echo "All experiments launched. Attach with: tmux attach -t benchmarks_${experiments[0]}_${WEIGHT_MODE}"
