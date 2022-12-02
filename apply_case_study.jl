# import packages
include("load_packages.jl")
include("functions/call_benchmark.jl")

# specify the weeks of training
all_training_days = 7:7:35
# specify the name of the datasource
datasource = "casestudy"
# specify the capacity of the warehouses in APPs
capacity_app = [12000,12000]
# specify the capacity of the warehouses in categorybrands
capacity_cb = [2550,6900]
# specify alpha for the heuristic
sig = 0.01
# specify whether to simulate on aggregated or unaggregated data
aggregated = true
# specify whether aall articles over the three month are allocated to warehouses or whether
# to use a rolling horizon of all_training_days + 14 training_days
rollingarticles = true
# specify whether to simulate the parcel output based on warehouse data
sim_ware = false
# specify whether to simulate the parcel output based on the order data
sim_order = true

# choose Optimisations and Heuristics to evaluate in the benchmark
start = DataFrame(
    QMKO  = [0], # quadratic-multiple knapsack heuristic with CPLEX as solver
    QMK   = [0], # quadratic-multiple knapsack heuristic with SBB as solver
    CHIM  = [1], # main chi-square heuristic without local search
    CHI   = [1], # chi-square heuristic + local search based on the QMK objective function
    KL    = [0], # K-LINK heuristic by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2020.08.024
    KLQ   = [0], # K-LINK optimisation with SBB by Zhang, W.-H. Lin, M. Huang and X. Hu (2021) https://doi.org/10.1016/j.ejor.2020.08.024
    GO    = [0], # greedy orders heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
    GP    = [0], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
    GS    = [1], # greedy seeds heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
    BS    = [1], # bestselling heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
    OPT   = [0]  # optimisation model to determine the optimal solution with CPLEX
)

# theme and colours
theme(
    :default,
    size=(800, 450),
    xlabel="date",
    margin=5mm,
    lw=2.0,
)
sc3 = palette(["steelblue3","darkgoldenrod1","gray10"])
scwh2 = palette(["firebrick","cyan4"])
scwh4 = palette(["salmon","firebrick","paleturquoise3","cyan4"])

# functions
function dict_eval(dict, key)
     get.(Ref(dict), key, 0)
end

function warehouse_weight(W,brands)
    weight = zeros(Int64,axes(W,1),axes(W,2))
    for i in axes(W,1)
        for j in axes(W,2)
            if W[i,j] == true
                weight[i,j] = round(Int64,brands[i,:articles])
            end
        end
    end
    weight = sum(weight,dims=1)
    return weight
end

function reorders(X,W,dict_algorithm,categorybrands,algorithm)
    ReorderCategorybrands = sum(X[:,:,dict_algorithm[algorithm]]) - sum(X[:,:,dict_algorithm[algorithm]] .* W[:,:])
    ReorderApps = sum(sum(X[:,:,dict_algorithm[algorithm]], dims = 2) .* categorybrands[:,:articles])
    ReorderApps = round(ReorderApps - sum(sum(X[:,:,dict_algorithm[algorithm]] .* W[:,:], dims = 2) .* categorybrands[:,:articles]))
    X[:,:,dict_algorithm[algorithm]] = copy(W[:,:])
    return ReorderCategorybrands, ReorderApps, X
end

function reorders(X,W)
    ReorderApps = sum(X) - sum(X .* W)
    X = copy(W)
    return ReorderApps, X
end

# import the warehouse data 
warehouse   = CSV.read("casestudy_data/warehouse_detail.csv", DataFrame)

# import the order data
orders     = CSV.read("casestudy_data/orders.csv", DataFrame)

# save the available yearweeks
yearweeks = unique(warehouse.yearweek) == unique(orders.yearweek) ? unique(warehouse.yearweek) : error("Warehouse and Orders do not match for yearweeks!")
alldates = unique(warehouse.date)

no_cat = length(unique(warehouse[:,:category]))
no_brand = length(unique(warehouse[:,:brand]))
no_size = length(unique(warehouse[:,:size]))

# prepare the available articles per category and brand
articlecategorybrands = groupby(warehouse, [:article,:category,:brand,:size,:date])
articlecategorybrands = combine(articlecategorybrands, :warehouse_A => mean => :warehouse_A, :warehouse_B => mean => :warehouse_B)
articlecategorybrands[!,:warehouse_A] = [articlecategorybrands[x,:warehouse_A] > 0 ? 1 : 0 for x in axes(articlecategorybrands,1)]
articlecategorybrands[!,:warehouse_B] = [articlecategorybrands[x,:warehouse_B] > 0 ? 1 : 0 for x in axes(articlecategorybrands,1)]

# plot the unique articles in each warehouse per week
plot_wh_aggregated = combine(groupby(articlecategorybrands,[:date]),nrow => :articles, :warehouse_A => sum => :warehouse_A, :warehouse_B => sum => :warehouse_B)
display(plot(plot_wh_aggregated.date, [plot_wh_aggregated.articles, plot_wh_aggregated.warehouse_A, plot_wh_aggregated.warehouse_B], labels=["articles" "warehouse_A" "warehouse_B"], ylabel = "articles", xlabel ="date"))

# aggregate the articles to categorybrands
categorybrands = combine(groupby(articlecategorybrands, [:category,:brand,:date]), nrow => :articles, :warehouse_A => sum => :warehouse_A, :warehouse_B => sum => :warehouse_B)
categorybrands = combine(groupby(categorybrands, [:category,:brand]), :articles => mean => :articles, :warehouse_A => mean => :warehouse_A, :warehouse_B => mean => :warehouse_B)
sort!(categorybrands, :articles, rev=true)

# create some dicts for the evaluation
dict_article_category = Dict(articlecategorybrands[x,:article] => articlecategorybrands[x,:category] for x in axes(articlecategorybrands,1))
dict_article_categorybrand = Dict(articlecategorybrands[x,:article] => (articlecategorybrands[x,:category],articlecategorybrands[x,:brand]) for x in axes(articlecategorybrands,1))
dict_categorybrand_id = Dict((categorybrands[x,:category],categorybrands[x,:brand]) => x for x in axes(categorybrands,1))
dict_categorybrand_weight = Dict((categorybrands[x,:category],categorybrands[x,:brand]) => categorybrands[x,:articles] for x in axes(categorybrands,1))

# estimate real split-delivery ratio from casestudy based on the orders alone
order_products = combine(groupby(orders, [:order,:warehouse_id]), nrow => :products)
order_parcels = combine(groupby(order_products, [:order]), nrow => :parcels)
casestudy_real_split_day = (sum(order_parcels[:,:parcels]) - nrow(order_parcels))/length(groupby(warehouse, [:date]))
casestudy_real_ratio_day = sum(order_parcels[:,:parcels])/nrow(order_parcels)-1

print("\n\nOverall Results from the Supplied Order Data",
    "\n   Average Split Deliveries per Day: ", casestudy_real_split_day,
    "\n   Average Split Ratio per Day:      ", casestudy_real_ratio_day)

# estimate real split_delivery ratio based on dispatch simulation
function real_splits(warehouse,orders,capacity_app)
    print("\n\nStarting Simulation of Order Dispatching based on Warehouse Data.")
    simparcels = DataFrame(data=String[], articles=Int64[], orders_test=Int64[], day=Date[],
            min_dispatch_A=Int64[], min_dispatch_B=Int64[], max_dispatch_A=Int64[], max_dispatch_B=Int64[],
            weight_app_A=Int64[], weight_app_B=Int64[], mode=String[], parcel_test=Int64[], reorderapps = Float64[],
    )
    unique_articles = unique(warehouse[:,:article])
    dict_article_id = Dict(unique_articles[x] => x for x in axes(unique_articles,1))
    orders_day = groupby(orders, [:date])
    warehouse_day = groupby(warehouse, [:date])
    X = zeros(Bool, length(unique_articles), 2)
    split_deliveries = 0
    changes = 0
    for testday in 1:length(warehouse_day)
        unique_test_ordernumbers = unique(orders_day[testday][:,:order])
        unique_article_warehouse = nrow(unique(select(orders_day[testday], [:article, :warehouse_id])))
        dict_test_orders_id = Dict(unique_test_ordernumbers[x] => x for x in axes(unique_test_ordernumbers,1))
        test_orders = sparse(dict_eval(dict_test_orders_id,orders_day[testday].order),dict_eval(dict_article_id,orders_day[testday].article),true)
        test_orders = [test_orders spzeros(Bool,length(unique_test_ordernumbers),length(unique_articles)-size(test_orders,2))]

        W = zeros(Bool,length(unique_articles),2)
        warehouse_local = warehouse_day[testday]
        for wh = 1:2
            for row = 1:nrow(warehouse_local)
                if wh == 1
                    W[dict_article_id[warehouse_local[row,:article]],wh] = warehouse_local[row,:warehouse_A]
                else
                    W[dict_article_id[warehouse_local[row,:article]],wh] = warehouse_local[row,:warehouse_B]
                end
            end
        end
        articlecapacity = vec(sum(W,dims=1))
        combination = COMBINEWAREHOUSES(capacity_app)
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(test_orders, W, articlecapacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(test_orders, W, articlecapacity, combination, false)
        ReorderApps, X = reorders(X,W)
        split_deliveries += parcels_benchmark
        changes += ReorderApps
        print("\n Day: ",orders_day[testday][1,:date]," / Splitdeliveries: ",parcels_benchmark, " / Reorders: ", ReorderApps,
            " / A: ", split_bench_max[1]," to ",split_bench_min[1]," / B: ", split_bench_min[2]," to ",split_bench_max[2],
            " / W: ", sum(W), " / T: ", unique_article_warehouse)
        push!(simparcels, (
            data=datasource, 
            articles=length(unique_article_warehouse), 
            orders_test=length(unique_test_ordernumbers), 
            day=orders_day[testday][1,:date], 
            min_dispatch_A=split_bench_max[1],
            min_dispatch_B=split_bench_min[2],
            max_dispatch_A=split_bench_min[1],
            max_dispatch_B=split_bench_max[2],
            weight_app_A=sum(W,dims=1)[1], 
            weight_app_B=sum(W,dims=1)[2], 
            mode="warehousesimulation", 
            parcel_test=parcels_benchmark, 
            reorderapps = ReorderApps,
        ))
    end
    casestudy_sim_split_day = split_deliveries / length(warehouse_day)
    casestuy_sim_ratio_day = sum(split_deliveries + nrow(order_parcels)) / nrow(order_parcels) -1
    casestudy_sim_reorders_day = changes/ length(warehouse_day)
    print("\n\nOverall Results from Simulation based on Warehousedata ",
    "\n   Average Split Deliveries per Day: ", casestudy_sim_split_day,
    "\n   Average Split Ratio per Day:      ", casestuy_sim_ratio_day,
    "\n   Average Reorders per Day:         ", casestudy_sim_reorders_day, )
    CSV.write("casestudy_results/warehousesim.csv",simparcels)
    return simparcels
end
sim_ware == true ? simparcels = real_splits(warehouse,orders,capacity_app) : nothing

if sim_ware == true
    display(@df simparcels plot(:day, 
        :parcel_test, ylabel="split deliveries", 
        label="Shipment Simulation",
        color_palette = sc3,
        title = "Split Deliveries: Shipment Simulation"
    ))

    display(@df simparcels plot(:day, 
    :parcel_test./:orders_test, ylabel="split ratio", 
    label="Shipment Simulation",
    color_palette = sc3,
    title = "Split Ratio: Shipment Simulation"
    ))

    display(@df simparcels plot(:day, 
        :reorderapps, ylabel="daily app movement", 
        label="Shipment Simulation",
        color_palette = sc3,
        title = "Transshipments: Shipment Simulation"
    ))

    display(@df simparcels plot(:day, 
        [:min_dispatch_A./:orders_test,:max_dispatch_A./:orders_test,:min_dispatch_B./:orders_test,:max_dispatch_B./:orders_test],
        ylabel="% Orders dispatched from Warehouse", 
        label=["Minimum A" "Maximum A" "Minimum B" "Maximum B"],
        color_palette = scwh4,
        title = "Flexibility: Shipment Simulation"
    ))

    display(@df simparcels plot(:day, 
        [:weight_app_A,:weight_app_B],
        ylabel="Warehouse Capacity used in APPs", 
        label=["Warehouse A" "Warehouse B"],
        color_palette = scwh2,
        title = "Capacity: Shipment Simulation"
    ))
end

# prepare data for weekly analysis based on optimisation
function benchmark_splits(start,aggregated,rollingarticles,all_training_days,alldates,categorybrands,orders,datasource,capacity_app,capacity_cb)
    print("\n\nStarting Optimisation of Warehouse Allocation based on Order Data.")
    benchmark = DataFrame(
        data=String[],
        skus=Int64[],
        orders_train=Int64[],
        orders_test=Int64[],
        date=Date[],
        traindays=Int64[],
        capacity_A=Int64[],
        capacity_B=Int64[],
        min_dispatch_A=Int64[],
        min_dispatch_B=Int64[],
        max_dispatch_A=Int64[],
        max_dispatch_B=Int64[],
        weight_app_A=Int64[],
        weight_app_B=Int64[],
        diff=Float64[],
        buffer=Float64[],
        mode=String[],
        parcel_train=Int64[],
        parcel_test=Int64[],
        reordercategorybrands = Int64[],
        reorderapps = Float64[],
        duration=Float64[],
        cap_used=Int64[],
        local_search=Int64[],
        gap=Float64[],
    )

    aggregated == true ? capacity = capacity_cb : capacity = capacity_app
    
    for training_days in all_training_days
        X = zeros(Bool, nrow(categorybrands), 2, size(start,2))

        for daynumber in training_days+1:length(alldates)

            print("\n\n Trainingdays: ", training_days)
            print("\n Testday: ", alldates[daynumber])
            train_orders_raw = filter(row -> row.date in alldates[daynumber-training_days:daynumber-1], orders)
            test_orders_raw = filter(row -> row.date == alldates[daynumber], orders)

            if rollingarticles == true
                unique_articles = unique(filter(row -> row.date in alldates[daynumber-training_days:min(daynumber+14,daynumber)], orders)[:,:article])
            else
                unique_articles = unique(warehouse[:,:article])
            end
            dict_article_id = Dict(unique_articles[x] => x for x in axes(unique_articles,1))

            unique_train_ordernumbers = unique(train_orders_raw[:,:order])
            unique_test_ordernumbers = unique(test_orders_raw[:,:order])

            dict_train_orders_id = Dict(unique_train_ordernumbers[x] => x for x in axes(unique_train_ordernumbers,1))
            dict_test_orders_id = Dict(unique_test_ordernumbers[x] => x for x in axes(unique_test_ordernumbers,1))

            if aggregated == true
                train_orders = sparse(dict_eval(dict_train_orders_id,train_orders_raw.order),dict_eval(dict_categorybrand_id,dict_eval(dict_article_categorybrand,train_orders_raw.article)),true)
                train_orders = [train_orders spzeros(Bool,length(unique_train_ordernumbers),nrow(categorybrands)-size(train_orders,2))]

                test_orders = sparse(dict_eval(dict_test_orders_id,test_orders_raw.order),dict_eval(dict_categorybrand_id,dict_eval(dict_article_categorybrand,test_orders_raw.article)),true)
                test_orders = [test_orders spzeros(Bool,length(unique_test_ordernumbers),nrow(categorybrands)-size(test_orders,2))]
            else
                train_orders = sparse(dict_eval(dict_train_orders_id,train_orders_raw.order),dict_eval(dict_article_id,train_orders_raw.article),true)
                train_orders = [train_orders spzeros(Bool,length(unique_train_ordernumbers),length(unique_articles)-size(train_orders,2))]

                test_orders = sparse(dict_eval(dict_test_orders_id,test_orders_raw.order),dict_eval(dict_article_id,test_orders_raw.article),true)
                test_orders = [test_orders spzeros(Bool,length(unique_test_ordernumbers),length(unique_articles)-size(test_orders,2))]
            end

            benchmark!(
                start,
                X,
                benchmark,
                categorybrands,
                datasource,
                training_days,
                alldates[daynumber],
                train_orders, 
                test_orders, 
                capacity, 
                1800, 
                false, 
                8,
                0.00,
                10000000,
                sig,
                false,
                50,
                100,
                10,
                2,
                false
            )
        end
    end
    CSV.write("casestudy_results/benchmark_$(sig)_$(minimum(all_training_days))_to_$(maximum(all_training_days))_$(capacity[1])_$(capacity[2]).csv",benchmark)
    return benchmark
end

warehouse = nothing
articlecategorybrands = nothing

sim_order == true ? simbench = benchmark_splits(start,aggregated,rollingarticles,all_training_days,alldates,categorybrands,orders,datasource,capacity_app,capacity_cb) : nothing

if sim_order == true
    display(@df filter(row -> row.mode in ["CHIM_0.01","BS","RND"], simbench) plot(:date, 
        :parcel_test ,groups=:mode,ylabel="split deliveries", 
        label=["Competing Heuristic" "Our Heuristic" "Random"],
        color_palette = sc3,
        title = "Split Deliveries: Compared Algorithms"
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01", "BS"], simbench) plot(:date, 
        :parcel_test ,groups=:mode,ylabel="split deliveries", 
        label=["Competing Heuristic" "Our Heuristic"],
        color_palette = sc3,
        title = "Split Deliveries: Compared Algorithms"
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01","BS","RND"], simbench) plot(:date, 
        :parcel_test./:orders_test ,groups=:mode,ylabel="split ratio", 
        label=["Competing Heuristic" "Our Heuristic" "Random"],
        color_palette = sc3,
        title = "Split Ratio: Compared Algorithms"
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01","BS"], simbench) plot(:date, 
        :parcel_test./:orders_test ,groups=:mode,ylabel="split ratio", 
        label=["Competing Heuristic" "Our Heuristic"],
        color_palette = sc3,
        title = "Split Ratio: Compared Algorithms"
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01","BS"], simbench) plot(:date, 
        :reorderapps ,groups=:mode,ylabel="daily app movement", 
        label=["Competing Heuristic" "Our Heuristic"],
        color_palette = sc3,
        title = "Transshipments: Compared Algorithms "
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01"], simbench) plot(:date, 
        [:min_dispatch_A./:orders_test,:max_dispatch_A./:orders_test,:min_dispatch_B./:orders_test,:max_dispatch_B./:orders_test],
        ylabel="% Orders dispatched from Warehouse", 
        label=["Minimum A" "Maximum A" "Minimum B" "Maximum B"],
        color_palette = scwh4,
        title = "Flexibility: Our Algorithm"
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01"], simbench) plot(:date, 
        [:weight_app_A,:weight_app_B],
        ylabel="Warehouse Capacity used in APPs", 
        label=["Warehouse A" "Warehouse B"],
        color_palette = scwh2,
        title = "Capacity: Our Algorithm"
    ))

    display(@df filter(row -> row.mode in ["BS"], simbench) plot(:date, 
        [:weight_app_A,:weight_app_B],
        ylabel="Warehouse Capacity used in APPs", 
        label=["Warehouse A" "Warehouse B"],
        color_palette = scwh2,
        title = "Capacity: Competing Algorithm"
    ))
end