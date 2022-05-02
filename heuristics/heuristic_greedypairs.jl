# greedy pairs heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function GREEDYPAIRS(Q::Matrix{Int64},
                     capacity::Array{Int64,1})
    if CHECKCAPACITY(Q::Array{Int64,2}, capacity::Array{Int64,1}) == 1
        ## Calculate the coapperance matrix
        capacity_left = copy(capacity)
        ## Sort the pairs of distinct SKUs by decreasing coappearance
        pairs = SORTPAIRS(Q::Array{Int64,2})
        X = Array{Int64,2}(undef,size(Q,2),size(capacity,1)) .= 0
        GREEDYPAIRSMAIN!(pairs::Array{Int64,2}, X::Array{Int64,2}, capacity_left::Array{Int64,1})
        ## While there is unallocated space in the DCs do
        X = FILLUP(X::Array{Int64,2}, Q::Array{Int64,2}, capacity_left::Array{Int64,1})
        return X::Array{Int64,2}
    end
end