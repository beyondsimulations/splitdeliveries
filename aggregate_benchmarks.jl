# Script to aggregate the results of all benchmarks into one file
## import packages
include("load_packages.jl")

using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics

## the following aggregates all results
experiments     = ["small_benchmark","medium_benchmark","large_benchmark"]
dependencies    = ["HD-HF","MD-HF","ID-HF","HD-SF","MD-SF","ID-SF"]

function load_data()
    frame = DataFrame[]
    for experiment in experiments
        for dependency in dependencies
            loadframe = CSV.read("results/$(experiment)_$dependency.csv", DataFrame)
            isempty(frame) ? frame = loadframe : frame = append!(frame,loadframe)
        end
    end
    return frame
end

frame = load_data()
CSV.write("results/overall_results.csv",frame)
