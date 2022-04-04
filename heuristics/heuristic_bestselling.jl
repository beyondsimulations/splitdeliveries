 # bestselling heuristic by A. Catalan and M. Fisher (2012) doi:10. 2139/ssrn.2166687
function BESTSELLING(trans::Array{Int64,2},
                     Q::Matrix{Int64},
                     capacity::Array{Int64,1})
    if CHECKCAPACITY(trans::Array{Int64,2},
                     capacity::Array{Int64,1}) == 1
        #  Start the heuristic
        capacity_left = copy(capacity::Array{Int64,1})
        B = BESTSELLING_B(size(trans,2),capacity_left)
        sales = SORTSALES(trans::Array{Int64,2})
        X = Array{Int64,2}(undef,size(trans,2),size(capacity,1)) .= 0
        # while there exists a DC d such that the overall capacity is smaller than B do
        # for such each DC d do
        BESTSELLINGSTART!(sales::Array{Int64,2},
                          capacity_left::Array{Int64,1},
                          B::Int64,
                          X::Array{Int64,2})
        if B == 0
            for j = 1:size(capacity_left,1)
                if capacity_left[j] > 0
                    capacity_left[j] -= 1
                end
            end
            X[sales[1,1],1] = 1
            ## Mark the SKU that is now allocated and the capacity used
            sales[1,4] = 1
            ## Start the assignment of the other SKUs
            GREEDYSEEDSTART!(sales::Array{Int64,2}, X::Array{Int64,2}, Q::Array{Int64,2}, capacity_left::Array{Int64,1})
            ## For each SKUs sorted by decreasing sales that has not been allocated yet
            GREEDYSEEDMAIN!(sales::Array{Int64,2}, X::Array{Int64,2}, Q::Array{Int64,2}, capacity_left::Array{Int64,1})
            ## While there is unallocated space in the DCs do
            X = FILLUP(X::Array{Int64,2}, Q::Array{Int64,2}, capacity_left::Array{Int64,1})
        else
            ## Assign the B top selling SKUs to each warehouse with capacity > B
            BESTSELLINGTOP!(sales::Array{Int64,2}, capacity_left::Array{Int64,1}, B::Int64, X::Array{Int64,2})
            GREEDYSEEDMAIN!(sales::Array{Int64,2}, X::Array{Int64,2}, Q::Array{Int64,2}, capacity_left::Array{Int64,1})
            X = FILLUP(X::Array{Int64,2}, Q::Array{Int64,2}, capacity_left::Array{Int64,1})
        end
        return X::Array{Int64,2}
    end
end