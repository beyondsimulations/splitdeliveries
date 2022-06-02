function CHISQUAREHEUR(trans::SparseMatrixCSC{Bool, Int64},
                       capacity::Array{Int64,1},
                       Q::Array{Int64,2},
                       sig::Float64,
                       localsearch::Bool)
    if CHECKCAPACITY(Q,capacity) == 1
        # Number of transactions
        J = size(trans,1)

        # Number of products
        I = size(Q,2)

        # Chi Square Test
        ## We start the chi-square heuristic by performing chi-square tests 
        ## for all unique SKU-combinations. For each undirected SKU-pair we 
        ## evaluate a chi-square test and either accept or reject the null 
        ## hypothesis of independence. Our goal is a matrix dep with the 
        ## “net-benefit” of all unique positive dependent SKU-combinations 
        ## and a matrix nor with all "independent contributions" of the unique 
        ## SKU-combinations. More details in our article.
        sum_cond_sku = Int.(vec(sum(trans,dims=1)))
        nor    = zeros(Float64,I,I)
        dep    = zeros(Float64,I,I)
        HYOPTHESISCHI!(nor,dep,Q,I,J,sig,sum_cond_sku)

        ## Create the sku-warehouse allocation matrix
        X = zeros(Bool,I,size(capacity,1))

        ## Determine the sum of the coappearances for each SKU on the base of
        ## the matrices dep and nor
        sum_dep = vec(sum(dep,dims = 1))
        sum_nor = vec(sum(nor,dims = 1))

        ## Calculate the weight of each warehouse. It shows us the density of 
        ## the independent coappearances in each warehouse if we were to allocate
        ## all SKUs simply according to the highest independent coappearances.
        weight = WHWEIGHT(capacity,sum_nor)

        ## Copy the capacity to have an array that keeps the left-over capacity
        cap_left = copy(capacity)

        ## Create an array that saves the current dependencies for all unallocated
        ## SKUs to all warehouses and a vector that holds all unallocated skus
        state_dep = zeros(Float64,size(X,1),length(cap_left))
        state_nor = zeros(Float64,size(X,1),length(cap_left))
        allocated = zeros(Bool,size(X,1))

        ## In the heuristic we compare the split-delivery minimising potential 
        ## of SKU allocations to maximise dependent coappearances (SKUs with 
        ## positive dependencies are stored in the same warehouse) with the 
        ## maximisation of independent coappearances (SKUs with high independent 
        ## coappearances are stored in warehouses with high capacities).
        ## Start the allocation:
        ## If all SKUs are allocated once, continue after the while loop, else
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
        i,k  = SELECTIK(sum_dep,sum_nor,weight,cap_left,X,dep,allocated)
        ALLOCATEONE!(X,dep,nor,sum_dep,sum_nor,state_dep,state_nor,cap_left,allocated,i,k)

        ## Check for all unallocated SKUs whether they have positive 
        ## dependencies to the SKUs in the warehouse k the last SKU was allocated 
        ## to. If so, check whether the dependencies are expected to dominate the 
        ## independent coapperances. If yes, allocate the corresponding SKUs to 
        ## the warehouse k.
        ADDDEPENDENT!(X,cap_left,k,dep,nor,sum_dep,sum_nor,state_dep,state_nor,allocated)
        FILLLAST!(X,cap_left,allocated)
        end
        if sum(cap_left) > 0
            ## First, remove every so far assigned dependent SKU-pair from the coappearance 
            ## matrix dep to prevent the allocation bias described in our article. 
            #Qs = REMOVEALLOC(X::Array{Bool,2},
            #                cap_left::Array{Int64,1},
            #                nor::Array{Float64,2},
            #                dep::Array{Float64,2})
            
            ## Afterwards allocate the SKUs with the highest potential allocation value to 
            ## each warehouse with leftover storage space until it is full. If no SKU 
            ## with coappearances is found and there is still storage space left, terminate the 
            ## algorithm as further allocations pose no benefit. 
            FILLUP!(X,Q,cap_left)
        end
    end
    if localsearch == true
        LOCALSEARCHCHI!(X,Q)
    end
    X = convert(Matrix{Int64},X)
    ## return the resulting allocation matrix
    return X::Array{Int64,2}
end







    
        



