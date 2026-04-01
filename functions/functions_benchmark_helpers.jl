function RUN_LS!(
    benchmark,
    W,
    trans_train,
    trans_test,
    Q_ls,
    capacity,
    combination,
    max_ls,
    mode_name,
    dependency,
    skus,
    wareh,
    diff,
    buffer,
    wmode,
    benchnr,
    orders,
    train_test,
    time_Q_ls,
    time_alg,
)
    W_ls = copy(W)
    sleep(0.01)
    GC.gc()
    time_ls =
        time_alg +
        time_Q_ls +
        @elapsed ls_count = LOCALSEARCHCHI!(
            trans_train, W_ls, Q_ls, capacity, false, 0, max_ls
        )
    parcels_test = PARCELSSEND(trans_test, W_ls, capacity, combination)
    parcels_train = PARCELSSEND(trans_train, W_ls, capacity, combination)
    flex = FLEXIBILITY(trans_test, W_ls)
    ls_mode = "$(mode_name)_LS"
    print(
        "\n    ",
        ls_mode,
        ": parcels test data: ",
        parcels_test,
        " / parcels training data: ",
        parcels_train,
        " / flex: ",
        flex,
        " / time: ",
        round(time_ls; digits = 3),
        " / local search: ",
        ls_count,
        " / warehouse: ",
        sum(W_ls; dims = 1),
    )
    push!(
        benchmark,
        (
            dependency = dependency,
            skus = skus,
            wareh = wareh,
            diff = diff,
            buffer = buffer,
            weight_mode = wmode,
            mode = ls_mode,
            benchiter = benchnr,
            orders = orders,
            train_test = train_test,
            parcel_train = parcels_train,
            parcel_test = parcels_test,
            flexibility = flex,
            duration = time_ls,
            cap_used = sum(W_ls),
            local_search = ls_count,
            gap = 0,
        ),
    )
end
