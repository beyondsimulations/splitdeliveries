# greedy seeds heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function GREEDYSEEDS(trans::SparseMatrixCSC{Bool, Int64},
                     Q::Matrix{Int64},
                     capacity::Array{Int64,1})
    if CHECKCAPACITY(Q,capacity) == 1
        ## Calculate the coapperance matrix
        capacity_left = copy(capacity) .- 1
        ## Sort the SKUs by decreasing sales
        sales = SORTSALES(trans)
        ## Assign the top (best selling) SKU to the largest DC
        X = Array{Bool,2}(undef,size(trans,2),size(capacity,1)) .= 0
        X[sales[1,1],1] = 1
        ## Mark the SKU that is now allocated and the capacity used
        sales[1,4] = 1
        ## Start the assignment of the other SKUs
        GREEDYSEEDSTART!(sales,X,Q,capacity_left)
        ## For each SKUs sorted by decreasing sales that has not been allocated yet
        GREEDYSEEDMAIN!(sales,X,Q,capacity_left)
        ## While there is unallocated space in the DCs do
        FILLUP!(X,Q,capacity_left)
        ## Return the solution
        X = convert(Matrix{Int64},X)
        return X::Array{Int64,2}
    end
end



    