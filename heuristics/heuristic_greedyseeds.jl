# greedy seeds heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function GREEDYSEEDS(
    trans::SparseMatrixCSC{Bool,Int64}, 
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real},
    )
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)
    # Create Coapperance Matrix
    Q = COAPPEARENCE(trans,sku_weight)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(capacity,sku_weight)
        ## Calculate the coapperance matrix
        capacity_left::Vector{Float64} = copy(capacity)
        ## Sort the SKUs by decreasing sales
        sales = SORTSALES(trans)
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
        # Check whether all SKUs are allocated
        if any(y->y < 1,sum(X,dims=2))
            "\n Error: Not all SKUs are allocated."
        end
        ## Return the solution
        return X
    end
end



    