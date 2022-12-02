# greedy pairs heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function GREEDYPAIRS(trans::SparseMatrixCSC{Bool,Int64},
                     capacity::Array{Int64,1},
                     sku_weight::Vector{<:Real})
    # Sort the warehouses by decreasing capacity
    #capacity = sort(capacity, rev=true)
    # Create Coapperance Matrix
    Q = COAPPEARENCE(trans,sku_weight)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(capacity,sku_weight)
        ## Calculate the coapperance matrix
        capacity_left::Vector{Float64} = copy(capacity)
        ## Sort the pairs of distinct SKUs by decreasing coappearance
        pairs = SORTPAIRS(Q)
        X = zeros(Bool,size(Q,2),size(capacity,1))
        GREEDYPAIRSMAIN!(pairs,X,capacity_left,sku_weight)
        ## While there is unallocated space in the DCs do
        FILLUP!(X,Q,capacity_left,sku_weight)
        # Check whether all SKUs are allocated
        if any(y->y < 1,sum(X,dims=2))
            "\n Error: Not all SKUs are allocated."
        end
        return X
    end
end