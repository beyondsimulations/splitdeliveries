using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics
using LaTeXStrings

experiments     = ["2_100skus","3_1000skus","4_10000skus"]
dependencies    = ["IND","MD","HD"]
datasets        = ["a_parcels_sent","b_duration"]

function load_data()
    frame = DataFrame[]
    for experiment in experiments
        for dataset in datasets
            for dependency in dependencies
                loadframe = CSV.read("results/$(experiment)_$(dataset)_$dependency.csv", DataFrame)
                loadcapacity = CSV.read("results/$(experiment)_$(dataset)_$dependency.csv", DataFrame)
                insertcols!(loadframe, 1, :experiment .=> "$(experiment)")
                insertcols!(loadframe, 2, :dataset .=> "$(dataset)")
                insertcols!(loadframe, 3, :dependency .=> "$dependency")
                insertcols!(loadframe, 6, :capacity_base .=> 0)
                insertcols!(loadframe, 7, :capvariation .=> 0)
                insertcols!(loadframe, 9, :buffer_rel .=> 0.0)
                loadframe[:,:capacity_base] = loadframe[:,:capacity] - loadframe[:,:buffer]
                for row = 1:nrow(loadframe)
                    loadframe[row,:buffer] == 0.0 ? nothing : loadframe[row,:buffer_rel] = loadframe[row,:buffer]/loadframe[row,:capacity_base]
                end
                allowmissing!(loadframe)
                base  = loadframe[1,:wareh]
                base_ver = 0
                for i = 1:nrow(loadframe)
                    if loadframe[i,:wareh] == base
                        base_ver += 1
                        loadframe[i,:capvariation] = base_ver
                    else
                        base_ver = 1
                        base = loadframe[i,:wareh]
                        loadframe[i,:capvariation] = base_ver
                    end
                end
                for i = 10:ncol(loadframe)
                    for j = 1:nrow(loadframe)
                        loadframe[j,i] == 0.0 ? loadframe[j,i] = missing : nothing
                    end
                end
                loadframe = stack(loadframe,10:ncol(loadframe))
                isempty(frame) ? frame = loadframe : frame = append!(frame,loadframe)
            end
        end
    end
    return frame
end

frame = load_data()
CSV.write("results/aggregated.csv",frame)

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