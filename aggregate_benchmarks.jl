# Script to aggregate the results of all benchmarks into one file
## import packages
include("load_packages.jl")

using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics

## the following aggregates all results
experiments = ["100_benchmark", "1000_benchmark", "10000_benchmark", "100000_benchmark"]
dependencies = ["HD-VF", "MD-VF", "ID-VF", "HD-SF", "MD-SF", "ID-SF"]
weight_suffixes = ["", "_frequency", "_random"]

function load_data(folder::String)
    frame = DataFrame[]
    for experiment in experiments
        for dependency in dependencies
            for suffix in weight_suffixes
                file = "$folder/$(experiment)_$dependency$suffix.csv"
                isfile(file) || continue
                loadframe = CSV.read(file, DataFrame; stringtype = String)
                isempty(frame) ? frame = loadframe : frame = append!(frame, loadframe)
            end
        end
    end
    return frame
end

is_chi(mode) = startswith(mode, "CHI")

# The main benchmark run predates the independence gate, so its CHI rows are
# replaced with the gated CHI results from the dedicated re-run. All other
# heuristics and the RND baseline stem from the main run.
frame = load_data("results")
frame_gated = load_data("results_chigate")

frame_merged = vcat(frame[.!is_chi.(frame.mode), :], frame_gated[is_chi.(frame_gated.mode), :])
CSV.write("results/overall_results.csv", frame_merged)

# Ungated CHI rows from the main run, used for the gate ablation.
CSV.write("results/overall_results_chiungated.csv", frame[is_chi.(frame.mode), :])
