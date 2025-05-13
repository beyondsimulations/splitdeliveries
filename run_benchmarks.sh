#!/bin/bash

# ./run_benchmarks.sh

# Array of dependencies
dependencies=(
    "ID-SF"
    "ID-VF"
    #"MD-SF"
    #"MD-VF"
    #"HD-SF"
    #"HD-VF"
)

experiment="k"

# Create a new tmux session (if not already in one)
session_name="benchmarks_${experiment}"
tmux new-session -d -s $session_name

# For each dependency, create a new window and run the benchmark
for i in "${!dependencies[@]}"; do
    dep="${dependencies[$i]}"

    # Create a new window
    tmux new-window -t $session_name -n "$dep"

    # Create the temporary Julia file for this dependency
    cat > "run_${dep}.jl" << EOL
include("load_packages.jl")
dependencies=["${dep}"]
experiment = "${experiment}"
include("main_benchmark_settings.jl")
EOL

    # Send the command to run Julia with this configuration
    tmux send-keys -t $session_name:$dep "julia run_${dep}.jl" C-m
done

# Attach to the tmux session if not already in tmux
if [ -z "$TMUX" ]; then
    tmux attach-session -t $session_name
fi
