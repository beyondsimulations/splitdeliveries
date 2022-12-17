 # bestselling heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function BESTSELLING(
    trans::SparseMatrixCSC{Bool,Int64}, 
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real}
    )
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)
    # Create Coapperance Matrix
    Q = COAPPEARENCE(trans,sku_weight)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(capacity,sku_weight)
        #  Start the heuristic
        capacity_left::Vector{Float64} = copy(capacity)
        sales = SORTSALES(trans)
        X = zeros(Bool,size(trans,2),size(capacity,1))
        # apply the Greedy Seeds Heuristic if B = 0
        if sum(capacity) == size(trans,2)
            ## Assign the top (best selling) SKU to the largest DC
            X = zeros(Bool,size(trans,2),size(capacity,1))
            X[sales[1,1],1] = 1
            ## Mark the SKU that is now allocated and the capacity used
            sales[1,4] = 1
            capacity_left[1] -= sku_weight[sales[1,1]]
            ## Start the assignment of the other SKUs
            GREEDYSEEDSTART!(sales,X,Q,capacity_left,sku_weight)
            ## For each SKUs sorted by decreasing sales that has not been allocated yet
            GREEDYSEEDMAIN!(sales,X,Q,capacity_left,sku_weight)
            ## While there is unallocated space in the DCs do
            FILLUP!(X,Q,capacity_left,sku_weight)
            ## Return the solution
        else
            B = BESTSELLING_B(capacity_left,sku_weight)
            # while there exists a DC d such that the overall capacity is smaller than B do
            # for such each DC d do
            BESTSELLINGSTART!(sales,capacity_left,B,X,sku_weight)
            ## Assign the B top selling SKUs to each warehouse with capacity > B
            BESTSELLINGTOP!(sales,capacity_left,capacity,B,X,sku_weight)
            GREEDYSEEDMAIN!(sales,X,Q,capacity_left,sku_weight)
            FILLUP!(X,Q,capacity_left,sku_weight)
        end
        # Check whether all SKUs are allocated
        if any(y->y < 1,sum(X,dims=2))
            "\n Error: Not all SKUs are allocated."
        end
        return X
    end
end