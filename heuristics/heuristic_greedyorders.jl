function GREEDYORDERS(
    trans::SparseMatrixCSC{Bool,Int64}, 
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real}
    )
    # Sort the warehouses by decreasing capacity
    #capacity = sort(capacity, rev=true)
    # Compute the Coapperance Matrix
    Q = COAPPEARENCE(trans,sku_weight)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(capacity,sku_weight)
        # Count the number of SKUs per order
        SKU_in_order = vec(sum(trans, dims=2))
        # Sort the orders by decreasing SKUs per order
        order_sort = sortperm(SKU_in_order,rev=true)
        # Create the allocation matrix
        X = zeros(Bool,size(Q,2),size(capacity,1))
        # Create a vector saving the leftover capacity of each warehouse
        capacity_left::Vector{Float64} = copy(capacity)
        # Start the main part of the greedy orders heuristic
        GREEDYORDERSMAIN!(order_sort,trans,capacity_left,X,sku_weight)
        ## While there is unallocated space in the DCs do
        FILLUP!(X,Q,capacity_left,sku_weight)
        # Check whether all SKUs are allocated
        if any(y->y < 1,sum(X,dims=2))
            "\n Error: Not all SKUs are allocated."
        end
        return X
    end
end