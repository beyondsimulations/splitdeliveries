
    GC.gc()

    ## import packages
    include("load_packages.jl")

    ## generate transactions
    print("\nstarting generation of transactions.")
    @time trans = RANDOMTRANS(40000,1000000,1000,0.00,0.01,0.00,0.50,0.10,0.10)

let
    ## import packages
    include("load_packages.jl")
    
    # state the capacites
    capacity = [15000,15000,10000]

    ## call the heuristic chi
    print("starting the heuristic chi")
    @time W = CHISQUAREHEUR(trans,capacity,0.01,true,true)
    print("heuristic finished")

    ## benchmark the parcels sent out
    combination = COMBINEWAREHOUSES(capacity)
    parcels = PARCELSSEND(trans,W,capacity,combination)
    print("\nParcels chi: ",parcels,"\n")

    ## call the heuristic
    print("starting the heuristic bs")
    @time W = BESTSELLING(trans,capacity)
    print("heuristic finished")

    ## benchmark the parcels sent out
    combination = COMBINEWAREHOUSES(capacity)
    parcels = PARCELSSEND(trans,W,capacity,combination)
    print("\nParcels gs: ",parcels,"\n")
end

sleep(1)
GC.gc()
