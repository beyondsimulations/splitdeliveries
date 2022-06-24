using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics
using LaTeXStrings

experiments = ["2_100skus","3_1000skus","4_10000skus"]
dependencies = ["IND","MD","HD"]

function load_data(experiment, dependencies)
    parcels = DataFrame[]
    duration = DataFrame[]
    for dependency in dependencies
        loadframe_parcels = CSV.read("results/$(experiment)_a_parcels_sent_$dependency.csv", DataFrame)
        loadframe_duration = CSV.read("results/$(experiment)_b_duration_$dependency.csv", DataFrame)
        insertcols!(loadframe_parcels, 1, :dependency .=> "$dependency")
        insertcols!(loadframe_duration, 1, :dependency .=> "$dependency")
        loadframe_duration = select!(loadframe_duration, Not(:RND))
        isempty(parcels) ? parcels = loadframe_parcels : parcels = append!(parcels,loadframe_parcels)
        isempty(duration) ? duration = loadframe_duration : duration = append!(duration,loadframe_duration)
    end
    insertcols!(parcels, 4, :capacity_base .=> 0)
    insertcols!(parcels, 6, :buffer_rel .=> 0.0)
    insertcols!(duration, 4, :capacity_base .=> 0)
    insertcols!(duration, 6, :buffer_rel .=> 0.0)
    parcels[:,:capacity_base] = parcels[:,:capacity] - parcels[:,:buffer]
    duration[:,:capacity_base] = duration[:,:capacity] - duration[:,:buffer]
    for row = 1:nrow(parcels)
        parcels[row,:buffer] == 0.0 ? nothing : parcels[row,:buffer_rel] = parcels[row,:buffer]/parcels[row,:capacity_base]
        duration[row,:buffer] == 0.0 ? nothing : duration[row,:buffer_rel] = duration[row,:buffer]/duration[row,:capacity_base]
    end
    allowmissing!(parcels)
    allowmissing!(duration)
    for i = 7:ncol(duration)
        for j = 1:nrow(duration)
            duration[j,i] == 0.0 ? duration[j,i] = missing : nothing
            parcels[j,i] == 0.0 ? parcels[j,i] = missing : nothing
        end
    end
    parcels = stack(parcels,7:ncol(parcels))
    duration = stack(duration,7:ncol(duration))
    return parcels, duration
end


for experiment in experiments
    parcels, duration = load_data(experiment, dependencies)

    duration_skus = groupby(dropmissing(duration),[:dependency,:capacity_base,:variable])
    duration_skus = combine(duration_skus, :value => mean => :value)
    for x = 1:length(dependencies)
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

end