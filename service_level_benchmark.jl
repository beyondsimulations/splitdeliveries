## import packages
include("load_packages.jl")

# Choose the benchmark which should be evaluated
## Benchmarks used in our article:
### b1_20skus
### b2_100skus
### b3_1000skus
### b4_10000skus
## Dependencies used in our article:
### EF
### ID
### MD
### HD
    experiment = "s1_1000skus"
    dependency = "ID"

#  Specify the number of orders and the ratio between test
## and training data for the generated transactional data sets
    orders     = 200000
    train_test = 0.50

# load the data that specifies the dependencies
    include("dependency/$dependency.jl")

# Set the number of cpu cores your computer has at its disposal
    cpu_cores  = 8
    ren_lock = ReentrantLock()


# Parameters for CHISQUARE
## sig: significance level alpha for the chi-square tests
    service_level = vcat(1.0:-0.1:0.2,0.1:-0.01:0.02,0.01:-0.001:0.002,0.001:-0.0001:0.0002,0.0001:-0.00001:0.00002,0.00001:-0.000001:0.000002)

# Parameters for RANDOM
## iterations: number of different random allocations for the comparison
    iterations = 100

# Initialise the basic problem by loading the respective capacity constellations
## capacity_benchmark: capacity matrix with column = capacity and row = constellation
    capacity_benchmark  = readdlm("capacity/capacity_$experiment.csv", ';', Int64)
    skus_benchmark      = vec(readdlm("capacity/skus_$experiment.csv", ';', Int64))
    
# Start the benchmark
    benchmark = DataFrame(dependency = String[], skus = Int64[], service_level = Float64[], wareh= Int64[], diff = Float64[], 
                            buffer = Float64[], heuristic = String[],
                            parcel_train = Int64[], parcel_test = Int64[], duration = Float64[])
    for trials = 1:10
        # Create the transactional data sets
        time = @elapsed trans = RANDOMTRANS(1000,orders,skus_in_order,sku_frequency,
                                            ceil(Int64,max(1000/50,10)),
                                            min_dependence,max_dependence,
                                            group_link,ind_chance,one_direction,
                                            multi_relatio)
        print("\ntransactions generated after ", round(time,digits = 3)," seconds.")

        #  Split the data into training and test data
        if train_test > 0.00
            cut = round(Int64,size(trans,1) * train_test)
            trans_train = trans[1:cut,:]
            trans_test  = trans[(cut+1):size(trans,1),:]
        else
            trans_train = trans_test = trans
        end

        for a = 1:size(capacity_benchmark,1)        
            ## Load the capacity of each individual run
            ### Note: It has to be sorted starting with the largest capacity
            capacity = Array{Int64,1}(undef,count(x -> x > 0, capacity_benchmark[a,:]))
            for b = 1:size(capacity,1)
                capacity[b] = capacity_benchmark[a,b]
            end
            print("\ncapacity constellation: ",a," of ",size(capacity_benchmark,1),
                "\ncapacity: ",capacity)

            for sig in service_level
                ## Create all possible capacity combinations for the parcel Benchmark
                combination = COMBINEWAREHOUSES(capacity)

                print("\nBenchmark for significance level of ",sig,".")
                ## Start chi square heuristic without local search
                sleep(0.1)
                GC.gc()
                time_benchmark = @elapsed W = CHISQUAREHEUR(trans_train,capacity,sig,false,false)
                parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                print("\n   CHIM: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))
                
                push!(benchmark, (dependency = dependency, skus = skus_benchmark[a], service_level = sig, wareh = length(capacity), 
                                diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                                heuristic = "CHIM", parcel_train = parcels_train, parcel_test = parcels_benchmark, 
                                duration = time_benchmark))

                ## Start chi square heuristic with local search
                sleep(0.1)
                GC.gc()
                time_benchmark = @elapsed W = CHISQUAREHEUR(trans_train,capacity,sig,true,false)
                parcels_benchmark = PARCELSSEND(trans_test, W, capacity, combination)
                parcels_train = PARCELSSEND(trans_train, W, capacity, combination)
                print("\n    CHI: parcels test data: ", parcels_benchmark, 
                        " / parcels training data: ", parcels_train,  
                        " / time: ",round(time_benchmark, digits = 3))

                push!(benchmark, (dependency = dependency, skus = skus_benchmark[a], service_level = sig, wareh = length(capacity), 
                diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                heuristic = "CHI", parcel_train = parcels_train, parcel_test = parcels_benchmark, 
                duration = time_benchmark))

                ## Benchmark the random allocation of SKUs
                sleep(0.1)
                time_benchmark = @elapsed parcels_benchmark = RANDOMBENCH(trans_test,capacity,iterations,combination)
                parcels_train = RANDOMBENCH(trans_train,capacity,iterations,combination)
                print("\n    RND: parcels test data: ", parcels_benchmark,
                    " / parcels training data: ", parcels_train,"\n")

                push!(benchmark, (dependency = dependency, skus = skus_benchmark[a], service_level = sig, wareh = length(capacity), 
                diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), buffer = round((sum(capacity)/skus_benchmark[a])-1,digits = 2),
                heuristic = "RND", parcel_train = parcels_train, parcel_test = parcels_benchmark, 
                duration = time_benchmark))

                CSV.write("results/$(experiment)_service_level_$dependency.csv", benchmark)
            end
        end
    end                                             
print("\nbenchmark finished at ",now(),".")
    