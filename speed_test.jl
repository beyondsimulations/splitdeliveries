
    GC.gc()

    ## import packages
    include("load_packages.jl")

    ## generate transactions
    print("\nstarting generation of transactions.")
    @time trans = RANDOMTRANS(20000,1000000,1000,0.00,0.01,0.00,0.50,0.10,0.10)

let
    ## import packages
    include("load_packages.jl")
    
    # state the capacites
    capacity = [8000,6000,6000]

    ## call the heuristic chi
    print("starting the heuristic chi")
    @time W = CHISQUAREHEUR(trans,capacity,0.01,true,true)
    print("heuristic finished")

    ## call the heuristic qmk
    print("starting the heuristic qmk")
    @time W = MQKP(trans,capacity,3600,"SBB",true,8,0.000,100000000000,"QMK")
    print("heuristic finished")
    

    ## call the heuristic
    print("starting the heuristic gp")
    @time W = GREEDYPAIRS(trans,capacity)
    print("heuristic finished")
end

sleep(1)
GC.gc()
