#!/bin/bash

# ./run_benchmarks.sh
# Runs experiments 100 and 1000 in parallel, each in its own tmux session.
# Within each experiment, all 6 dependencies run in parallel (tmux windows).

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
experiments=("100" "1000")

for experiment in "${experiments[@]}"; do
    session_name="benchmarks_${experiment}"
    tmux new-session -d -s $session_name

    for i in "${!dependencies[@]}"; do
        dep="${dependencies[$i]}"

        # Create a new window
        tmux new-window -t $session_name -n "$dep"

        # Use experiment-specific filename to avoid conflicts between parallel runs
        cat > "run_${experiment}_${dep}.jl" << EOL
    include("load_packages.jl")
    dependencies=["${dep}"]
    experiment = "${experiment}"
    include("main_benchmark_settings.jl")
EOL

        tmux send-keys -t $session_name:$dep "julia run_${experiment}_${dep}.jl" C-m
    done

    echo "Started experiment=${experiment} in tmux session '${session_name}'"
done

echo "All experiments launched. Attach with: tmux attach -t benchmarks_100"
