# greedy pairs heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function GREEDYPAIRS(Q::Matrix{Int64},
                     capacity::Array{Int64,1})
    if CHECKCAPACITY(Q::Array{Int64,2}, capacity::Array{Int64,1}) == 1
        ## Calculate the coapperance matrix
        capacity_left = copy(capacity)
        ## Sort the pairs of distinct SKUs by decreasing coappearance
        pairs = SORTPAIRS(Q)
        X = Array{Bool,2}(undef,size(Q,2),size(capacity,1)) .= 0
        GREEDYPAIRSMAIN!(pairs,X,capacity_left)
        ## While there is unallocated space in the DCs do
        FILLUP!(X,Q,capacity_left)
        X = convert(Matrix{Int64},X)
        return X::Array{Int64,2}
    end
end