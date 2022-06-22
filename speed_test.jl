
    GC.gc()

    ## import packages
    include("load_packages.jl")

    ## generate transactions
    print("\nstarting generation of transactions.")
    @time trans = RANDOMTRANS(50000,10000000,1000,0.00,0.01,0.00,0.50,0.10,0.10)

let
    # state the capacites
    capacity = [17000,17000,16000]

    ## call the heuristic chi
    print("starting the heuristic chi")
    @time W = CHISQUAREHEUR(trans,capacity,0.01,true,true)
    print("heuristic finished")

    ## call the heuristic qmk
    print("starting the heuristic bs")
    @time W = BESTSELLING(trans,capacity)
    print("heuristic finished")
end

sleep(1)
GC.gc()
