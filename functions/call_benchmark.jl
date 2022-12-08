function benchmark!(
    start::DataFrame,
    X::Array{Bool,3},
    benchmark::DataFrame,
    categorybrands::DataFrame,
    datasource::String,
    training_days::Int64,
    date,
    trans_train, 
    trans_test, 
    capacity, 
    abort, 
    show_opt, 
    cpu_cores,
    allowed_gap,
    max_nodes,
    sig,
    chistatus,
    max_ls,
    trials,
    stagnant,
    strategy,
    klinkstatus,
    sku_weight,
    optimization_cycle::Bool,
    training_intervall,
    )

## Determine the possible warehouse combinations
    combination = COMBINEWAREHOUSES(capacity)

## Determine dict to allocate the necessary X
    dict_algorithm = Dict(names(start)[i] => i for i in axes(start,2))

## Start QMK heuristic to find the optimal solution with the solver CPLEX
if start[1,:QMKO] == 1
    if all(y->y == 1,sku_weight)
        sleep(0.01)
        GC.gc()
        if optimization_cycle == true
            time_benchmark = @elapsed W,gap_optimisation = MQKP(trans_train,capacity,sku_weight,abort,"CPLEX",show_opt,
                                                                cpu_cores,allowed_gap,max_nodes,"QMK")
        else
            W = X[:,:,dict_algorithm["QMKO"]]
            time_benchmark = 0.0
            gap_optimisation = 0.0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus, X = reorders(X,W,dict_algorithm,categorybrands,"QMKO")
        print("\n   QMKO: parcels test data: ", parcels_benchmark, 
            " / parcels training data: ", round(Int64, parcels_train/training_days), 
            " / time: ",round(time_benchmark,digits = 3),
            " / gap: ",round(gap_optimisation,6),
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )
        
        push!(benchmark, (data = datasource, 
                            skus = size(trans_train,2),
                            orders_train = size(trans_train,1),
                            orders_test = size(trans_test,1),
                            date=date,
                            traindays=training_days,
                            trainingintervall=training_intervall,
                            capacity_47=capacity[1],
                            capacity_50=capacity[2],
                            min_dispatch_47=split_bench_max[1],
                            min_dispatch_50=split_bench_min[2],
                            max_dispatch_47=split_bench_min[1],
                            max_dispatch_50=split_bench_max[2],
                            weight_app_47=weight[1],
                            weight_app_50=weight[2],
                            diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                            buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                            mode = "QMKO",
                            parcel_train = round(Int64, parcels_train/training_days), 
                            parcel_test = parcels_benchmark,
                            reordercategorybrands = ReorderCategorybrands,
                            reorderskus = reorderskus,
                            duration = time_benchmark,
                            cap_used = sum(W),
                            local_search = 0,
                            gap = gap_optimisation))
    end
end

## Start QMK heuristic with SBB as solver
if start[1,:QMK] == 1
    if all(y->y == 1,sku_weight)
        sleep(0.01)
        GC.gc()
        if optimization_cycle == true
            time_benchmark = @elapsed W,gap_optimisation = MQKP(trans_train,capacity,sku_weight,abort, "SBB",show_opt,
                                                                cpu_cores,allowed_gap,max_nodes,"QMK")
        else
            W = X[:,:,dict_algorithm["QMK"]]
            time_benchmark = 0.0
            gap_optimisation = 0.0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus, X = reorders(X,W,dict_algorithm,categorybrands,"QMK")
        print("\n    QMK: parcels week test data: ", parcels_benchmark, 
            " / parcels week training data: ", round(Int64,parcels_train/training_days),  
            " / time: ", round(time_benchmark, digits = 3),
            " / gap: ",round(gap_optimisation, digits = 6),
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )

        push!(benchmark, (data = datasource, 
                        skus = size(trans_train,2),
                        orders_train = size(trans_train,1),
                        orders_test = size(trans_test,1),
                        date=date,
                        traindays=training_days,
                        trainingintervall=training_intervall,
                        capacity_47=capacity[1],
                        capacity_50=capacity[2],
                        min_dispatch_47=split_bench_max[1],
                        min_dispatch_50=split_bench_min[2],
                        max_dispatch_47=split_bench_min[1],
                        max_dispatch_50=split_bench_max[2],
                        weight_app_47=weight[1],
                        weight_app_50=weight[2],
                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                        buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                        mode = "QMK",
                        parcel_train = round(Int64, parcels_train/training_days), 
                        parcel_test = parcels_benchmark,
                        reordercategorybrands = ReorderCategorybrands,
                        reorderskus = reorderskus,
                        duration = time_benchmark,
                        cap_used = sum(W),
                        local_search = 0,
                        gap = gap_optimisation))
    end
end

## Start chi square heuristic without local search
if start[1,:CHIM] == 1
        sleep(0.01)
        GC.gc()
        if optimization_cycle == true
            time_benchmark = @elapsed W, ls = CHISQUAREHEUR(trans_train,capacity,sig,0,sku_weight,chistatus)
        else
            W = X[:,:,dict_algorithm["CHIM"]]
            time_benchmark = 0.0
            ls = 0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"CHIM")
        print("\n   CHIM: parcels test data: ", parcels_benchmark, 
            " / parcels training data: ", round(Int64, parcels_train/training_days),  
            " / time: ",round(time_benchmark, digits = 3),
            " / local search: ", ls,
            " / sig: ", sig,
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )

        push!(benchmark, (
                    data = datasource, 
                    skus = size(trans_train,2),
                    orders_train = size(trans_train,1),
                    orders_test = size(trans_test,1),
                    date=date,
                    traindays=training_days,
                    trainingintervall=training_intervall,
                    capacity_47=capacity[1],
                    capacity_50=capacity[2],
                    min_dispatch_47=split_bench_max[1],
                    min_dispatch_50=split_bench_min[2],
                    max_dispatch_47=split_bench_min[1],
                    max_dispatch_50=split_bench_max[2],
                    weight_app_47=weight[1],
                    weight_app_50=weight[2],
                    diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                    buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                    mode = "CHIM_$sig",
                    parcel_train = round(Int64, parcels_train/training_days), 
                    parcel_test = parcels_benchmark,
                    reordercategorybrands = ReorderCategorybrands,
                    reorderskus = reorderskus,
                    duration = time_benchmark,
                    cap_used = sum(W),
                    local_search = ls,
                    gap = 0))
end

## Start chi square heuristic with local search
if start[1,:CHI] == 1
    if all(y->y==sku_weight[1],sku_weight)
        sleep(0.01)
        GC.gc()
        if optimization_cycle == true
            time_benchmark = @elapsed W,ls = CHISQUAREHEUR(trans_train,capacity,sig,max_ls,sku_weight,chistatus)
        else
            W = X[:,:,dict_algorithm["CHI"]]
            time_benchmark = 0.0
            ls = 0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"CHI")
        print("\n    CHI: parcels test data: ", parcels_benchmark, 
            " / parcels training data: ", round(Int64, parcels_train/training_days),  
            " / time: ",round(time_benchmark, digits = 3),
            " / local search: ", ls,
            " / sig: ", sig,
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )

        push!(benchmark, (data = datasource, 
                    skus = size(trans_train,2),
                    orders_train = size(trans_train,1),
                    orders_test = size(trans_test,1),
                    date=date,
                    traindays=training_days,
                    trainingintervall=training_intervall,
                    capacity_47=capacity[1],
                    capacity_50=capacity[2],
                    min_dispatch_47=split_bench_max[1],
                    min_dispatch_50=split_bench_min[2],
                    max_dispatch_47=split_bench_min[1],
                    max_dispatch_50=split_bench_max[2],
                    weight_app_47=weight[1],
                    weight_app_50=weight[2],
                    diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                    buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                    mode = "CHI_$sig",
                    parcel_train = round(Int64, parcels_train/training_days), 
                    parcel_test = parcels_benchmark,
                    reordercategorybrands = ReorderCategorybrands,
                    reorderskus = reorderskus,
                    duration = time_benchmark,
                    cap_used = sum(W),
                    local_search = ls,
                    gap = 0))
    end
end

## Start our reproduction of the  K-LINKS heuristic by
## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
if  start[1,:KL] == 1
    if all(y->y==1,sku_weight)
        sleep(0.01)
        GC.gc()
        if optimization_cycle == true
            time_benchmark = @elapsed W,ls = KLINKS(trans_train,capacity,trials,stagnant,strategy,klinkstatus)
        else
            W = X[:,:,dict_algorithm["KL"]]
            time_benchmark = 0.0
            ls = 0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"KL")
        print("\n     KL: parcels train data: ", parcels_benchmark, 
            " / parcels training data: ", round(Int64, parcels_train/training_days),  
            " / time: ", round(time_benchmark, digits = 3),
            " / local search: ", ls,
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )

        push!(benchmark, (data = datasource, 
                        skus = size(trans_train,2),
                        orders_train = size(trans_train,1),
                        orders_test = size(trans_test,1),
                        date=date,
                        traindays=training_days,
                        trainingintervall=training_intervall,
                        capacity_47=capacity[1],
                        capacity_50=capacity[2],
                        min_dispatch_47=split_bench_max[1],
                        min_dispatch_50=split_bench_min[2],
                        max_dispatch_47=split_bench_min[1],
                        max_dispatch_50=split_bench_max[2],
                        weight_app_47=weight[1],
                        weight_app_50=weight[2],
                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                        buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                        mode = "KL",
                        parcel_train = round(Int64, parcels_train/training_days), 
                        parcel_test = parcels_benchmark,
                        reordercategorybrands = ReorderCategorybrands,
                        reorderskus = reorderskus,
                        duration = time_benchmark,
                        cap_used = sum(W),
                        local_search = ls,
                        gap = 0))
    end
end

## Start our reproduction of the  K-LINKS optimization with SBB by
## [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
if  start[1,:KLQ] == 1
    if all(y->y==1,sku_weight)
        sleep(0.01)
        GC.gc()
        if optimization_cycle
            time_benchmark = @elapsed W, gap_optimisation = MQKP(trans_train,capacity,sku_weight,abort,"SBB",show_opt,
                                                                cpu_cores,allowed_gap,max_nodes,"QMK")
        else
            W = X[:,:,dict_algorithm["KLQ"]]
            time_benchmark = 0.0
            gap_optimisation = 0.0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"KLQ")
        print("\n    KLQ: parcels test data: ", parcels_benchmark, 
            " / parcels training data: ", round(Int64, parcels_train/training_days),  
            " / time: ", round(time_benchmark, digits = 3),
            " / gap: ", round(gap_optimisation, digits = 6),
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )

        push!(benchmark, (data = datasource, 
                        skus = size(trans_train,2),
                        orders_train = size(trans_train,1),
                        orders_test = size(trans_test,1),
                        date=date,
                        traindays=training_days,
                        trainingintervall=training_intervall,
                        capacity_47=capacity[1],
                        capacity_50=capacity[2],
                        min_dispatch_47=split_bench_max[1],
                        min_dispatch_50=split_bench_min[2],
                        max_dispatch_47=split_bench_min[1],
                        max_dispatch_50=split_bench_max[2],
                        weight_app_47=weight[1],
                        weight_app_50=weight[2],
                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                        buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                        mode = "KLQ",
                        parcel_train = round(Int64, parcels_train/training_days), 
                        parcel_test = parcels_benchmark, 
                        reordercategorybrands = ReorderCategorybrands,
                        reorderskus = reorderskus,
                        duration = time_benchmark,
                        cap_used = sum(W),
                        local_search = 0,
                        gap = gap_optimisation))
    end
end

## Start our reproduction of the greedy orders heuristic by
## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
if start[1,:GO] == 1
    sleep(0.01)
    GC.gc()
    if optimization_cycle == true
        time_benchmark = @elapsed W = GREEDYORDERS(trans_train,capacity,sku_weight)
    else
        W = X[:,:,dict_algorithm["GO"]]
        time_benchmark = 0.0
    end
    parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
    parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
    parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
    weight = warehouse_weight(W,categorybrands)
    ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"GO")
    print("\n     GO: parcels test data: ", parcels_benchmark, 
        " / parcels training data: ", round(Int64, parcels_train/training_days),  
        " / time: ",round(time_benchmark, digits = 3),
        " \n        ",
        " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
        " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
        " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
        " / ReorderCB: ", ReorderCategorybrands,
        " / ReorderSKUs: ", reorderskus,
        )


    push!(benchmark, (data = datasource, 
                    skus = size(trans_train,2),
                    orders_train = size(trans_train,1),
                    orders_test = size(trans_test,1),
                    date=date,
                    traindays=training_days,
                    trainingintervall=training_intervall,
                    capacity_47=capacity[1],
                    capacity_50=capacity[2],
                    min_dispatch_47=split_bench_max[1],
                    min_dispatch_50=split_bench_min[2],
                    max_dispatch_47=split_bench_min[1],
                    max_dispatch_50=split_bench_max[2],
                    weight_app_47=weight[1],
                    weight_app_50=weight[2],
                    diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                    buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                    mode = "GO",
                    parcel_train = round(Int64, parcels_train/training_days), 
                    parcel_test = parcels_benchmark,
                    reordercategorybrands = ReorderCategorybrands,
                    reorderskus = reorderskus,
                    duration = time_benchmark,
                    cap_used = sum(W),
                    local_search = 0,
                    gap = 0))
end

## Start our reproduction of the greedy pairs heuristic by
## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
if start[1,:GP] == 1
    sleep(0.01)
    GC.gc()
    if optimization_cycle == true
        time_benchmark = @elapsed W = GREEDYPAIRS(trans_train,capacity,sku_weight)
    else
        W = X[:,:,dict_algorithm["GP"]]
        time_benchmark = 0.0
    end
    parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
    parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
    parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
    weight = warehouse_weight(W,categorybrands)
    ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"GP")
    print("\n     GP: parcels test data: ", parcels_benchmark, 
        " / parcels training data: ", round(Int64, parcels_train/training_days),  
        " / time: ",round(time_benchmark, digits = 3),
        " \n        ",
        " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
        " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
        " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
        " / ReorderCB: ", ReorderCategorybrands,
        " / ReorderSKUs: ", reorderskus,
        )
    
    push!(benchmark, (data = datasource, 
                    skus = size(trans_train,2),
                    orders_train = size(trans_train,1),
                    orders_test = size(trans_test,1),
                    date=date,
                    traindays=training_days,
                    trainingintervall=training_intervall,
                    capacity_47=capacity[1],
                    capacity_50=capacity[2],
                    min_dispatch_47=split_bench_max[1],
                    min_dispatch_50=split_bench_min[2],
                    max_dispatch_47=split_bench_min[1],
                    max_dispatch_50=split_bench_max[2],
                    weight_app_47=weight[1],
                    weight_app_50=weight[2],
                    diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                    buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                    mode = "GP",
                    parcel_train = round(Int64, parcels_train/training_days), 
                    parcel_test = parcels_benchmark,
                    reordercategorybrands = ReorderCategorybrands,
                    reorderskus = reorderskus,
                    duration = time_benchmark,
                    cap_used = sum(W),
                    local_search = 0,
                    gap = 0))
end

## Start our reproduction of the greedy seeds heuristic by
## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
if start[1,:GS] == 1
    sleep(0.01)
    GC.gc()
    if optimization_cycle == true
        time_benchmark = @elapsed W = GREEDYSEEDS(trans_train,capacity,sku_weight)
    else
        W = X[:,:,dict_algorithm["GS"]]
        time_benchmark = 0.0
    end
    parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
    parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
    parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
    weight = warehouse_weight(W,categorybrands)
    ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"GS")
    print("\n     GS: parcels test data: ", parcels_benchmark, 
        " / parcels training data: ", round(Int64, parcels_train/training_days),  
        " / time: ",round(time_benchmark, digits = 3),
        " \n        ",
        " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
        " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
        " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
        " / ReorderCB: ", ReorderCategorybrands,
        " / ReorderSKUs: ", reorderskus,
        )

    push!(benchmark, (data = datasource, 
                skus = size(trans_train,2),
                orders_train = size(trans_train,1),
                orders_test = size(trans_test,1),
                date=date,
                traindays=training_days,
                trainingintervall=training_intervall,
                capacity_47=capacity[1],
                capacity_50=capacity[2],
                min_dispatch_47=split_bench_max[1],
                min_dispatch_50=split_bench_min[2],
                max_dispatch_47=split_bench_min[1],
                max_dispatch_50=split_bench_max[2],
                weight_app_47=weight[1],
                weight_app_50=weight[2],
                diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                mode = "GS",
                parcel_train = round(Int64, parcels_train/training_days), 
                parcel_test = parcels_benchmark,
                reordercategorybrands = ReorderCategorybrands,
                reorderskus = reorderskus,
                duration = time_benchmark,
                cap_used = sum(W),
                local_search = 0,
                gap = 0))
end

## Start our reproduction of the  bestselling heuristic by
## [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
if  start[1,:BS] == 1
    sleep(0.01)
    GC.gc()
    if optimization_cycle
        time_benchmark = @elapsed W = BESTSELLING(trans_train,capacity,sku_weight)
    else
        W = X[:,:,dict_algorithm["BS"]]
        time_benchmark = 0.0
    end
    parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
    parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
    parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
    weight = warehouse_weight(W,categorybrands)
    ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"BS")
    print("\n     BS: parcels test data: ", parcels_benchmark, 
        " / parcels training data: ", round(Int64, parcels_train/training_days),  
        " / time: ",round(time_benchmark, digits = 3),
        " \n        ",
        " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
        " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
        " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
        " / ReorderCB: ", ReorderCategorybrands,
        " / ReorderSKUs: ", reorderskus,
        )
    
    push!(benchmark, (data = datasource, 
                skus = size(trans_train,2),
                orders_train = size(trans_train,1),
                orders_test = size(trans_test,1),
                date=date,
                traindays=training_days,
                trainingintervall=training_intervall,
                capacity_47=capacity[1],
                capacity_50=capacity[2],
                min_dispatch_47=split_bench_max[1],
                min_dispatch_50=split_bench_min[2],
                max_dispatch_47=split_bench_min[1],
                max_dispatch_50=split_bench_max[2],
                weight_app_47=weight[1],
                weight_app_50=weight[2],
                diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                mode = "BS",
                parcel_train = round(Int64, parcels_train/training_days), 
                parcel_test = parcels_benchmark,
                reordercategorybrands = ReorderCategorybrands,
                reorderskus = reorderskus,
                duration = time_benchmark,
                cap_used = sum(W),
                local_search = 0,
                gap = 0))
end

## Start the search for optimal solution with the solver CPLEX
## Choose FULLOPTEQ if each SKUs can only be allocated once, else use
## FULLOPTUEQ if SKUs can be allocated multiple times
if start[1,:OPT] == 1
    if all(y->y==1,sku_weight)
        sleep(0.01)
        GC.gc()
        if optimization_cycle
            if sum(capacity) == size(trans,2)
                time_benchmark = @elapsed W,gap_optimisation,popt = FULLOPTEQ(trans_train,capacity,abort,show_opt,
                                                                                cpu_cores,allowed_gap,max_nodes)
            else
                time_benchmark = @elapsed W,gap_optimisation,popt = FULLOPTUEQ(trans_train,capacity,abort,show_opt,
                                                                                cpu_cores,allowed_gap,max_nodes)
            end
        else
            W = X[:,:,dict_algorithm["OPT"]]
            time_benchmark = 0.0
            gap_optimisation = 0.0
            popt = 0.0
        end
        parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
        parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
        parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
        weight = warehouse_weight(W,categorybrands)
        ReorderCategorybrands, reorderskus = reorders(X,W,dict_algorithm,categorybrands,"OPT")
        print("\n    OPT: parcels test data: ", parcels_benchmark, 
            " / parcels training data: ", round(Int64, parcels_train/training_days),  
            " / time: ",round(time_benchmark, digits = 3),
            " \n        ",
            " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
            " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
            " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
            " / ReorderCB: ", ReorderCategorybrands,
            " / ReorderSKUs: ", reorderskus,
            )
        
        push!(benchmark, (data = datasource, 
                        skus = size(trans_train,2),
                        orders_train = size(trans_train,1),
                        orders_test = size(trans_test,1),
                        date=date,
                        traindays=training_days,
                        trainingintervall=training_intervall,
                        capacity_47=capacity[1],
                        capacity_50=capacity[2],
                        min_dispatch_47=split_bench_max[1],
                        min_dispatch_50=split_bench_min[2],
                        max_dispatch_47=split_bench_min[1],
                        max_dispatch_50=split_bench_max[2],
                        weight_app_47=weight[1],
                        weight_app_50=weight[2],
                        diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                        buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                        mode = "OPT",
                        parcel_train = round(Int64, parcels_train/training_days), 
                        parcel_test = parcels_benchmark,
                        reordercategorybrands = ReorderCategorybrands,
                        reorderskus = reorderskus,
                        duration = time_benchmark,
                        cap_used = sum(W),
                        local_search = 0,
                        gap = gap_optimisation))
    end
end

## Benchmark the random allocation of SKUs
sleep(0.01)
parcels_benchmark_rand = 0.0
split_bench_max_rand = [0.0,0.0]
split_bench_min_rand = [0.0,0.0]
weight_rand = [0.0,0.0]
time_benchmark = 0.0
for cycle in 1:10
    time_benchmark += @elapsed W = RANDOMALLOCMULTI(trans_test,capacity,sku_weight)
    parcels_benchmark, split_bench_max = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, true)
    parcels_benchmark, split_bench_min = PARCELSSEND_WEIGHT(trans_test, W, capacity, combination, false)
    parcels_train, split_train = PARCELSSEND_WEIGHT(trans_train, W, capacity, combination, true)
    weight = warehouse_weight(W,categorybrands)
    parcels_benchmark_rand += parcels_benchmark * (1/10)
    split_bench_max_rand .+= [split_bench_max[1] .*(1/10), split_bench_max[2] .*(1/10)]
    split_bench_min_rand .+= [split_bench_min[1] .*(1/10), split_bench_min[2] .*(1/10)]
    weight_rand .= [weight[1] .*(1/10), weight[2] .*(1/10)]
end
print("\n    RND: parcels test data: ", parcels_benchmark,
    " / parcels training data: ", round(Int64, parcels_train/training_days),
    " / time: ", round(time_benchmark, digits = 3),
    " \n        ",
    " / 47: ", split_bench_max[1]," to ",split_bench_min[1],
    " / 50: ", split_bench_min[2]," to ",split_bench_max[2],
    " / APP: 47 - ", weight[1]," and 50 -  ",weight[2],
    )

push!(benchmark, (data = datasource, 
                skus = size(trans_train,2),
                orders_train = size(trans_train,1),
                orders_test = size(trans_test,1),
                date=date,
                traindays=training_days,
                trainingintervall=training_intervall,
                capacity_47=capacity[1],
                capacity_50=capacity[2],
                min_dispatch_47=round(Int64,split_bench_max_rand[1]),
                min_dispatch_50=round(Int64,split_bench_min_rand[2]),
                max_dispatch_47=round(Int64,split_bench_min_rand[1]),
                max_dispatch_50=round(Int64,split_bench_max_rand[2]),
                weight_app_47=round(Int64,weight[1]),
                weight_app_50=round(Int64,weight[2]),
                diff = round((maximum(capacity)-minimum(capacity))/sum(capacity), digits = 2), 
                buffer = round((sum(capacity)/size(trans_train,2))-1,digits = 2),
                mode = "RND",
                parcel_train = round(Int64, parcels_train/training_days), 
                parcel_test = parcels_benchmark,
                reordercategorybrands = 0,
                reorderskus = 0.0,
                duration = time_benchmark,
                cap_used = sum(W),
                local_search = 0,
                gap = 0))
end