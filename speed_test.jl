using TimerOutputs
tmr = TimerOutput()
let
    ## import packages
    include("load_packages.jl")

    # state the capacites
    capacity = [10000,10000,10000]

    ## generate transactions
    print("\nstarting generation of transactions.")
    trans = RANDOMTRANS(30000,1000000,1000,0.00,0.01,0.00,0.50,0.10,0.10)

    ## generate the coappearance matrix
    @time @timeit tmr "coappearance matrix" Q = COAPPEARENCE(trans)

    ## call the heuristic
    print("\nstarting the heuristic")
    @time @timeit tmr "chisquare heuristic" W = CHISQUAREHEUR(trans, capacity, Q, 0.01, true)
    print("\nheuristic finished")

    ## benchmark the parcels sent out
    combination = COMBINEWAREHOUSES(capacity)
    parcels = PARCELSSEND(trans,W,capacity,combination)
    show(tmr)
    print("\n Parcels: ",parcels,"\n")
end

sleep(10)
GC.gc()
