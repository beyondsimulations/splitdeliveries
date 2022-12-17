## import packages
include("load_packages.jl")

using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics
#using LaTeXStrings

experiments     = ["e1_100to1000skus","e2_1000to10000skus","e3_10000to50000skus"]
dependencies    = ["HD-VF","MD-VF","ID-VF"]
datasets        = ["benchmark"]

function load_data()
    frame = DataFrame[]
    for experiment in experiments
        for dataset in datasets
            for dependency in dependencies
                loadframe = CSV.read("results/Published Results/$(experiment)_$(dataset)_$dependency.csv", DataFrame)
                isempty(frame) ? frame = loadframe : frame = append!(frame,loadframe)
            end
        end
    end
    return frame
end

frame = load_data()
CSV.write("results/overall_results.csv",frame)

#=for experiment in experiments
    parcels, duration = load_data()

    duration_skus = groupby(dropmissing(duration),[:dependency,:capacity_base,:variable])
    duration_skus = combine(duration_skus, :value => mean => :value)
    for x = 1:eachindex(dependencies)
        dependency = dependencies[x]
        frame = filter(:dependency => n -> n == "$(dependency)", duration_skus)
        plot(frame.capacity_base, frame.value, group = frame.variable, xlabel = "SKUs", ylabel = "Computation time", legend = :topleft, ylims = (0,500))
        savefig("graphs/duration_$(experiment)_$(dependency).pdf")
    end

    duration_skus = groupby(dropmissing(duration),[:capacity_base,:variable])
    duration_skus = combine(duration_skus, :value => mean => :value)
    frame = duration_skus
    display(plot(frame.capacity_base, frame.value, group = frame.variable, xlabel = "SKUs", ylabel = "Computation time", legend = :topleft, ylims = (0,500)))
    savefig("graphs/duration_$(experiment).pdf")


    pwb = groupby(dropmissing(parcels),[:dependency,:wareh,:buffer_rel,:variable])
    pwb = combine(pwb, :value => mean => :value)
    for dependency in dependencies
        frame = filter(:dependency => n -> n == "$dependency", pwb)
        display(bar((frame.wareh,frame.buffer_rel), frame.value, group = frame.variable, title = "$dependency", xlabel = "SKUs", ylabel = "Parcels dispatched", legend = :topleft))
    end

end =#