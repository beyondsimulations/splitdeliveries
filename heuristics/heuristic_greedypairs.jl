# greedy pairs heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function GREEDYPAIRS(trans::SparseMatrixCSC{Bool,Int64},
                     capacity::Array{Int64,1})
    # Create Coapperance Matrix
    Q = COAPPEARENCE(trans)
    # Clean the principle diagonal
    CLEANPRINCIPLE!(Q)
    if CHECKCAPACITY(Q,capacity) == 1
        ## Calculate the coapperance matrix
        capacity_left = copy(capacity)
        ## Sort the pairs of distinct SKUs by decreasing coappearance
        pairs = SORTPAIRS(Q)
        X = zeros(Bool,size(Q,2),size(capacity,1))
        GREEDYPAIRSMAIN!(pairs,X,capacity_left)
        ## While there is unallocated space in the DCs do
        FILLUP!(X,Q,capacity_left)
        return X
    end
end