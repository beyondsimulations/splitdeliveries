 # bestselling heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function BESTSELLING(trans::SparseMatrixCSC{Bool, Int64},
                     capacity::Array{Int64,1})
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)
    # Create Coapperance Matrix
    Q = COAPPEARENCE(trans)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(Q,capacity) == 1
        #  Start the heuristic
        capacity_left = copy(capacity)
        B = BESTSELLING_B(size(trans,2),capacity_left)
        sales = SORTSALES(trans)
        X = zeros(Bool,size(trans,2),size(capacity,1))
        # while there exists a DC d such that the overall capacity is smaller than B do
        # for such each DC d do
        BESTSELLINGSTART!(sales,capacity_left,B,X)
        if B == 0
            for j in axes(capacity_left,1)
                if capacity_left[j] > 0
                    capacity_left[j] -= 1
                end
            end
            X[sales[1,1],1] = 1
            ## Mark the SKU that is now allocated and the capacity used
            sales[1,4] = 1
            ## Start the assignment of the other SKUs
            GREEDYSEEDSTART!(sales,X,Q,capacity_left)
            ## For each SKUs sorted by decreasing sales that has not been allocated yet
            GREEDYSEEDMAIN!(sales,X,Q,capacity_left)
            ## While there is unallocated space in the DCs do
            FILLUP!(X,Q,capacity_left)
        else
            ## Assign the B top selling SKUs to each warehouse with capacity > B
            BESTSELLINGTOP!(sales,capacity_left,B,X)
            GREEDYSEEDMAIN!(sales,X,Q,capacity_left)
            FILLUP!(X,Q,capacity_left)
        end
        return X
    end
end