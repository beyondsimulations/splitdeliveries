# Script to aggregate the results of all benchmarks into one file
## import packages
include("load_packages.jl")

using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics

## the following aggregates all results
experiments     = ["small"]
dependencies    = ["ID-SF","ID-MF","ID-HF","MD-SF","MD-MF","MD-HF","HD-SF","HD-MF","HD-HF"]
datasets        = ["benchmark"]

function load_data()
    frame = DataFrame[]
    for experiment in experiments
        for dataset in datasets
            for dependency in dependencies
                loadframe = CSV.read("results/$(experiment)_$(dataset)_$dependency.csv", DataFrame)
                isempty(frame) ? frame = loadframe : frame = append!(frame,loadframe)
            end
        end
    end
    return frame
end

frame = load_data()
CSV.write("results/overall_results.csv",frame)
