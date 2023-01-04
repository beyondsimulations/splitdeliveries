function CHISQUAREHEUR(trans::SparseMatrixCSC{Bool,Int64},
                       capacity::Vector{Int64},
                       sig::Float64,
                       max_ls::Int64,
                       sku_weight::Vector{<:Real},
                       log_results::Bool)
    # Number of transactions
    J = size(trans,1)
    # Number of products
    I = size(trans,2)
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)
    # Create Coapperance Matrix
    log_results == true ? print("\n  creating coappearance matrix.") : nothing
    Q = COAPPEARENCE(trans,sku_weight)
    # Start the main part
    if CHECKCAPACITY(capacity,sku_weight)
        # Create Vector with number of transactions containing each SKU
        ordered_skus = [Q[i,i] for i in axes(Q,1)]
        # Clean the principle diagonal
        CLEANPRINCIPLE!(Q)
        # Chi Square Test
        ## We start the chi-square heuristic by performing chi-square tests 
        ## for all unique SKU-combinations. For each undirected SKU-pair we 
        ## evaluate a chi-square test and either accept or reject the null 
        ## hypothesis of independence. Our goal is a matrix dep with the 
        ## “net-benefit” of all unique positive dependent SKU-combinations 
        ## and a matrix nor with all "independent contributions" of the unique 
        ## SKU-combinations. More details in our article.
        log_results == true ? print("\n  starting chi-square tests.") : nothing
        dep          = zeros(Float64,I,I)
        HYOPTHESISTEST!(dep,Q,I,J,sig,ordered_skus)
        ## Create the sku-warehouse allocation matrix
        X = zeros(Bool,I,size(capacity,1))
        ## Determine the sum of the coappearances for each SKU on the base of
        ## the matrices dep and nor
        sum_dep = vec(sum(dep,dims = 2))./sku_weight
        sum_nor = (vec(sum(Q,dims = 2)) .- vec(sum(dep,dims = 2)))./sku_weight
        ## Copy the capacity to have an array that keeps the left-over capacity
        cap_left::Vector{Float64} = copy(capacity)
        ## Create arrays that save the current dependencies for all unallocated
        ## SKUs to all warehouses and a vector that holds all unallocated skus
        state_dep = zeros(Float64,size(X,1),length(cap_left))
        state_nor = zeros(Float64,size(X,1),length(cap_left))
        allocated = zeros(Bool,size(X,1))
        ## Calculate the weight of each warehouse. It shows us the density of 
        ## the independent coappearances in each warehouse if we were to allocate
        ## all SKUs simply according to the highest independent coappearances.
        weight = WHWEIGHT(capacity,sum_nor,sku_weight)
        ## In the heuristic we compare the split-delivery minimising potential 
        ## of SKU allocations to maximise dependent coappearances (SKUs with 
        ## positive dependencies are stored in the same warehouse) with the 
        ## maximisation of independent coappearances (SKUs with high independent 
        ## coappearances are stored in warehouses with high capacities).
        ## Start the allocation:
        ## If all SKUs are allocated once, continue after the while loop, else
        log_results == true ? print("\n  starting unique allocation of each sku.") : nothing
        ## Allocate all SKUs without coappearances to the smallest warehouses
        ## ALLOCATENOCOAPP!(X,dep,Q,sum_dep,sum_nor,state_dep,state_nor,cap_left,allocated,sku_weight)
        ## Allocate the rest of the SKUs with coapperances
        while sum(X) < I
        ## select the SKU i with the highest coappearance not being allocated
        ## to the warehouses yet. In addition, select the first sorted 
        ## warehouse k which has unused capacity left.
        ## Then, check for each warehouse with free space whether SKU i has 
        ## significant dependencies to other already allocated SKUs.
        ## If dependencies exist check whether they are expected to dominate.
        ## If so, allocate SKU i to to the warehouse with the dependencies 
        ## to complement the positive dependent SKUs in the warehouse. 
        ## Otherwise allocate SKU i to warehouse k as the independent 
        ## coappearances are expected to dominate.
        ## If no dependencies to yet allocated SKUs exist, select the warehouse
        ## with the highest unused capacity. If the dependet coappearances are
        ## expected to to dominate later, we allocate SKU i to the warehouse 
        ## with the highest unused capacity. This maximises the likelihood to 
        ## allocate all significant SKUs from a SKU-cluster into one warehouse. 
        ## Otherwise we allocate it to warehouse k to maximise the independent 
        ## coappearances in the allocation.
        i,k  = SELECTIK(sum_dep,sum_nor,weight,cap_left,X,dep,allocated,sku_weight)
        ALLOCATEONE!(X,dep,Q,sum_dep,sum_nor,state_dep,state_nor,cap_left,allocated,sku_weight,i,k)
        ## Check for all unallocated SKUs whether they have positive 
        ## dependencies to the SKUs in the warehouse k the last SKU was allocated 
        ## to. If so, check whether the dependencies are expected to dominate the 
        ## independent coapperances. If yes, allocate the corresponding SKUs to 
        ## the warehouse k.
        ADDDEPENDENT!(X,Q,cap_left,k,dep,sum_dep,sum_nor,state_dep,state_nor,allocated,sku_weight)
        FILLLAST!(X,cap_left,allocated,sku_weight)
        end
        if sum(cap_left) > 0
            ## First, remove every so far assigned dependent SKU-pair from the coappearance 
            ## matrix dep to prevent the allocation bias described in our article.
            log_results == true ? print("\n  starting allocation to leftover space.") : nothing
            REMOVEALLOC!(X,Q,dep)
            ## Afterwards allocate the SKUs with the highest potential allocation value to 
            ## each warehouse with leftover storage space until it is full. If no SKU 
            ## with coappearances is found and there is still storage space left, terminate the 
            ## algorithm as further allocations pose no benefit. 
            FILLUP!(X,dep,cap_left,sku_weight)
        end
        ## Free up RAM for last stage
        state_dep = nothing
        state_nor = nothing
        dep = nothing
        GC.gc()
        ls = 0
        if max_ls > 0
            if all(y->y==sku_weight[1],sku_weight)
                log_results == true ? print("\n  starting local search.") : nothing
                ls = LOCALSEARCHCHI!(trans,X,Q,capacity,log_results,ls,max_ls)
            else
                print("\n\n!!! Sorry, local search is only available if all SKUs have the same capacity.")
                print("\n!!! The following results will thus only reflect the CHIM heuristic.\n")
            end
        end
        # Check whether all SKUs are allocated
        if any(y->y < 1,sum(X,dims=2))
            "\n Error: Not all SKUs are allocated."
        end
    end
    ## return the resulting allocation matrix
    return X,ls
end
