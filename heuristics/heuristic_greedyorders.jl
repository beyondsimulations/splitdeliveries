function GREEDYORDERS(trans::SparseMatrixCSC{Bool,Int64}, capacity::Array{Int64,1})
    # Sort the warehouses by decreasing capacity
    capacity = sort(capacity, rev=true)
    # Compute the Coapperance Matrix
    Q = COAPPEARENCE(trans)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(Q,capacity) == 1
        # Count the number of SKUs per order
        SKU_in_order = vec(sum(trans, dims=2))
        # Sort the orders by decreasing SKUs per order
        order_sort = sortperm(SKU_in_order,rev=true)
        # Create the allocation matrix
        X = zeros(Bool,size(Q,2),size(capacity,1))
        # Create a vector saving the leftover capacity of each warehouse
        capacity_left = copy(capacity)
        # Start the main part of the greedy orders heuristic
        GREEDYORDERSMAIN!(order_sort,trans,capacity_left,X)
        ## While there is unallocated space in the DCs do
        FILLUP!(X,Q,capacity_left)
        return X
    end
end