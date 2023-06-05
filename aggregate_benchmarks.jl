## import packages
include("load_packages.jl")

using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics


## the following aggregates all results
experiments     = ["r1_100skus","r2_1000skus","r3_10000skus","r4_100000skus"]
dependencies    = ["HD-VF","MD-VF","ID-VF","HD-SF","MD-SF","ID-SF"]
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