# import packages
include("load_packages.jl")
include("functions/call_benchmark.jl")

# specify the weeks of training
all_training_days = 7:7
# specify the name of the datasource
datasource = "casestudy"
# specify the capacity of the warehouses in categorybrands
capacity = [400000,400000]
#capacity = [2550,6900]
# specify alpha for the heuristic
sig = 0.01
# specify whether to simulate binary capacity assignment or with group weights
binary = false
# specify whether aall articles over the three month are allocated to warehouses or whether
# to use a rolling horizon of all_training_days + 14 training_days
rollingarticles = true
# specify whether to simulate the parcel output based on warehouse data
sim_ware = true
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
    GP    = [1], # greedy pairs heuristic by A. Catalan and M. Fisher (2012) https://doi.org/10.2139/ssrn.2166687
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
    weight = zeros(Float64,axes(W,1),axes(W,2))
    for i in axes(W,1)
        for j in axes(W,2)
            if W[i,j] == true
                weight[i,j] = brands[i,:articles]
            end
        end
    end
    weight = floor.(sum(weight,dims=1))
    return weight
end

function reorders(X,W,dict_algorithm,categorybrands,algorithm)
    ReorderCategorybrands = sum(X[:,:,dict_algorithm[algorithm]]) - sum(X[:,:,dict_algorithm[algorithm]] .* W[:,:])
    ReorderSKUs = sum(sum(X[:,:,dict_algorithm[algorithm]], dims = 2) .* categorybrands[:,:articles])
    ReorderSKUs = round(ReorderSKUs - sum(sum(X[:,:,dict_algorithm[algorithm]] .* W[:,:], dims = 2) .* categorybrands[:,:articles]))
    X[:,:,dict_algorithm[algorithm]] = copy(W[:,:])
    return ReorderCategorybrands, ReorderSKUs, X
end

function reorders(X,W)
    ReorderSKUs = sum(X) - sum(X .* W)
    X = copy(W)
    return ReorderSKUs, X
end

# import the warehouse data 
warehouse = CSV.read("casestudy_data/warehouse_sku.csv", DataFrame)

# import the order data
orders = CSV.read("casestudy_data/orders.csv", DataFrame)

# import the order data
orders_sku = CSV.read("casestudy_data/orders_sku.csv", DataFrame)

# save the available yearweeks
yearweeks = unique(warehouse.yearweek) == unique(orders.yearweek) ? unique(warehouse.yearweek) : error("Warehouse and Orders do not match for yearweeks!")
alldates = unique(warehouse.date)

no_cat = length(unique(warehouse[:,:category]))
no_brand = length(unique(warehouse[:,:brand]))
no_size = length(unique(warehouse[:,:sku]))

# prepare the available articles per category and brand
articlecategorybrands = groupby(warehouse, [:article,:category,:brand,:sku,:date])
articlecategorybrands = combine(articlecategorybrands, :warehouse_47 => mean => :warehouse_47, :warehouse_50 => mean => :warehouse_50)
articlecategorybrands[!,:warehouse_47] = [articlecategorybrands[x,:warehouse_47] > 0 ? 1 : 0 for x in axes(articlecategorybrands,1)]
articlecategorybrands[!,:warehouse_50] = [articlecategorybrands[x,:warehouse_50] > 0 ? 1 : 0 for x in axes(articlecategorybrands,1)]

# plot the unique articles in each warehouse per week
plot_wh_aggregated = combine(groupby(articlecategorybrands,[:date]),nrow => :articles, :warehouse_47 => sum => :warehouse_47, :warehouse_50 => sum => :warehouse_50)
display(plot(plot_wh_aggregated.date, [plot_wh_aggregated.articles, plot_wh_aggregated.warehouse_47, plot_wh_aggregated.warehouse_50], labels=["articles" "warehouse_47" "warehouse_50"], ylabel = "articles", xlabel ="date"))

# aggregate the articles to categorybrands
categorybrands = combine(groupby(articlecategorybrands, [:category,:brand,:date]), nrow => :articles, :warehouse_47 => sum => :warehouse_47, :warehouse_50 => sum => :warehouse_50)
categorybrands = combine(groupby(categorybrands, [:category,:brand]), :articles => mean => :articles, :warehouse_47 => mean => :warehouse_47, :warehouse_50 => mean => :warehouse_50)
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
function real_splits(warehouse,orders,capacity)
    print("\n\nStarting Simulation of Order Dispatching based on Warehouse Data.")
    simparcels = DataFrame(data=String[], articles=Int64[], skus=Int64[], orders_test=Int64[], day=Date[],
            min_dispatch_47=Int64[], min_dispatch_50=Int64[], max_dispatch_47=Int64[], max_dispatch_50=Int64[],
            weight_app_47=Int64[], weight_app_50=Int64[], mode=String[], parcel_test=Int64[], parcel_real=Int64[], reorderskus = Float64[],
    )
    unique_articles = unique(warehouse[:,:article])
    unique_skus = unique(warehouse[:,:sku])
    dict_article_id = Dict(unique_articles[x] => x for x in axes(unique_articles,1))
    dict_sku_id = Dict(unique_skus[x] => x for x in axes(unique_skus,1))
    orders_day = groupby(orders_sku, [:date])
    warehouse_day = groupby(warehouse, [:date])
    X = zeros(Bool, length(unique_articles), 2)
    split_deliveries = 0
    changes = 0
    for testday in 1:length(warehouse_day)
        unique_test_ordernumbers = unique(orders_day[testday][:,:order])
        unique_article_warehouse = nrow(unique(select(orders_day[testday], [:article, :warehouse_id])))
        unique_sku_warehouse = nrow(unique(select(orders_day[testday], [:sku, :warehouse_id])))
        dict_test_orders_id = Dict(unique_test_ordernumbers[x] => x for x in axes(unique_test_ordernumbers,1))
        test_orders = sparse(dict_eval(dict_test_orders_id,orders_day[testday].order),dict_eval(dict_sku_id,orders_day[testday].sku),true)
        test_orders = [test_orders spzeros(Bool,length(unique_test_ordernumbers),length(unique_skus)-size(test_orders,2))]

        order_prod = combine(groupby(orders_day[testday], [:order,:warehouse_id]), nrow => :products)
        order_parcel = combine(groupby(order_prod, [:order]), nrow => :parcels)
        casestudy_real_split = (sum(order_parcel[:,:parcels]) - nrow(order_parcel))

        W = zeros(Bool,length(unique_skus),2)
        warehouse_local = warehouse_day[testday]
        for wh = 1:2
            for row = 1:nrow(warehouse_local)
                if wh == 1
                    W[dict_article_id[warehouse_local[row,:sku]],wh] = warehouse_local[row,:warehouse_47]
                else
                    W[dict_article_id[warehouse_local[row,:sku]],wh] = warehouse_local[row,:warehouse_50]
                end
            end
        end
        articlecapacity = vec(sum(W,dims=1))
        combination = COMBINEWAREHOUSES(capacity)
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(test_orders, W, articlecapacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(test_orders, W, articlecapacity, combination, false)
        ReorderSKUs, X = reorders(X,W)
        split_deliveries += parcels_benchmark
        changes += ReorderSKUs
        print("\n Day: ",orders_day[testday][1,:date]," / Splitdeliveries: ",parcels_benchmark, " / Reorders: ", ReorderSKUs,
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1]," / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / W: ", sum(W), " / T: ", unique_article_warehouse)
        push!(simparcels, (
            data=datasource, 
            articles=length(unique_article_warehouse),
            skus=length(unique_sku_warehouse),
            orders_test=length(unique_test_ordernumbers), 
            day=orders_day[testday][1,:date], 
            min_dispatch_47=split_bench_max[1],
            min_dispatch_50=split_bench_min[2],
            max_dispatch_47=split_bench_min[1],
            max_dispatch_50=split_bench_max[2],
            weight_app_47=sum(W,dims=1)[1], 
            weight_app_50=sum(W,dims=1)[2], 
            mode="warehousesimulation", 
            parcel_test=parcels_benchmark,
            parcel_real=casestudy_real_split,
            reorderskus = ReorderSKUs,
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
sim_ware == true ? simparcels = real_splits(warehouse,orders,capacity) : nothing

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
        :reorderskus, ylabel="daily app movement", 
        label="Shipment Simulation",
        color_palette = sc3,
        title = "Transshipments: Shipment Simulation"
    ))

    display(@df simparcels plot(:day, 
        [:min_dispatch_47./:orders_test,:max_dispatch_47./:orders_test,:min_dispatch_50./:orders_test,:max_dispatch_50./:orders_test],
        ylabel="% Orders dispatched from Warehouse", 
        label=["Minimum 47" "Maximum 47" "Minimum 50" "Maximum 50"],
        color_palette = scwh4,
        title = "Flexibility: Shipment Simulation"
    ))

    display(@df simparcels plot(:day, 
        [:weight_app_47,:weight_app_50],
        ylabel="Warehouse Capacity used in APPs", 
        label=["Warehouse 47" "Warehouse 50"],
        color_palette = scwh2,
        title = "Capacity: Shipment Simulation"
    ))
end

# prepare data for weekly analysis based on optimisation
function benchmark_splits(start,binary,rollingarticles,all_training_days,alldates,categorybrands,orders,datasource,capacity)
    print("\n\nStarting Optimisation of Warehouse Allocation based on Order Data.")
    benchmark = DataFrame(
        data=String[],
        skus=Int64[],
        orders_train=Int64[],
        orders_test=Int64[],
        date=Date[],
        traindays=Int64[],
        capacity_47=Int64[],
        capacity_50=Int64[],
        min_dispatch_47=Int64[],
        min_dispatch_50=Int64[],
        max_dispatch_47=Int64[],
        max_dispatch_50=Int64[],
        weight_app_47=Int64[],
        weight_app_50=Int64[],
        diff=Float64[],
        buffer=Float64[],
        mode=String[],
        parcel_train=Int64[],
        parcel_test=Int64[],
        reordercategorybrands = Int64[],
        reorderskus = Float64[],
        duration=Float64[],
        cap_used=Int64[],
        local_search=Int64[],
        gap=Float64[],
    )

    ## Determine the SKU weight
    if binary == true
        sku_weight = zeros(Int64,nrow(categorybrands)) .= 1
    else
        sku_weight = categorybrands[:,:articles]
    end
    
    for training_days in all_training_days
        X = zeros(Bool, nrow(categorybrands), 2, size(start,2))

        for daynumber in maximum(all_training_days)+1:length(alldates)

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

            train_orders = sparse(dict_eval(dict_train_orders_id,train_orders_raw.order),dict_eval(dict_categorybrand_id,dict_eval(dict_article_categorybrand,train_orders_raw.article)),true)
            train_orders = [train_orders spzeros(Bool,length(unique_train_ordernumbers),nrow(categorybrands)-size(train_orders,2))]

            test_orders = sparse(dict_eval(dict_test_orders_id,test_orders_raw.order),dict_eval(dict_categorybrand_id,dict_eval(dict_article_categorybrand,test_orders_raw.article)),true)
            test_orders = [test_orders spzeros(Bool,length(unique_test_ordernumbers),nrow(categorybrands)-size(test_orders,2))]

            #train_orders = sparse(dict_eval(dict_train_orders_id,train_orders_raw.order),dict_eval(dict_article_id,train_orders_raw.article),true)
            #train_orders = [train_orders spzeros(Bool,length(unique_train_ordernumbers),length(unique_articles)-size(train_orders,2))]

            #test_orders = sparse(dict_eval(dict_test_orders_id,test_orders_raw.order),dict_eval(dict_article_id,test_orders_raw.article),true)
            #test_orders = [test_orders spzeros(Bool,length(unique_test_ordernumbers),length(unique_articles)-size(test_orders,2))]

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
                false,
                sku_weight
            )
        end
    end
    CSV.write("casestudy_results/benchmark_$(sig)_$(minimum(all_training_days))_to_$(maximum(all_training_days))_$(capacity[1])_$(capacity[2]).csv",benchmark)
    return benchmark
end

warehouse = nothing
articlecategorybrands = nothing

sim_order == true ? simbench = benchmark_splits(start,binary,rollingarticles,all_training_days,alldates,categorybrands,orders,datasource,capacity) : nothing

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
        :reorderskus ,groups=:mode,ylabel="daily app movement", 
        label=["Competing Heuristic" "Our Heuristic"],
        color_palette = sc3,
        title = "Transshipments: Compared Algorithms "
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01"], simbench) plot(:date, 
        [:min_dispatch_47./:orders_test,:max_dispatch_47./:orders_test,:min_dispatch_50./:orders_test,:max_dispatch_50./:orders_test],
        ylabel="% Orders dispatched from Warehouse", 
        label=["Minimum 47" "Maximum 47" "Minimum 50" "Maximum 50"],
        color_palette = scwh4,
        title = "Flexibility: Our Algorithm"
    ))

    display(@df filter(row -> row.mode in ["CHIM_0.01"], simbench) plot(:date, 
        [:weight_app_47,:weight_app_50],
        ylabel="Warehouse Capacity used in APPs", 
        label=["Warehouse 47" "Warehouse 50"],
        color_palette = scwh2,
        title = "Capacity: Our Algorithm"
    ))

    display(@df filter(row -> row.mode in ["BS"], simbench) plot(:date, 
        [:weight_app_47,:weight_app_50],
        ylabel="Warehouse Capacity used in APPs", 
        label=["Warehouse 47" "Warehouse 50"],
        color_palette = scwh2,
        title = "Capacity: Competing Algorithm"
    ))
end