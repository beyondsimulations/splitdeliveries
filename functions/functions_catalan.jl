## function to sort the SKUs after the number of sales
function SORTSALES(trans::SparseMatrixCSC{Float64, Int64})
    sales = Array{Int64,2}(undef,size(trans,2),4) .= 0
    sales[:,3] = sum(trans,dims=1)
    sales[:,4] .= 0
    for i = 1:size(trans,2)
        sales[i,1] = i
    end
    sales = sortslices(sales,dims=1,by=x->x[3],rev=true)
    for i = 1:size(trans,2)
        sales[i,2] = i
    end
    # 1 column: index of SKU in allocation matrix
    # 2 column: index of SKU in sorted sales
    # 3 column: number of total sold SKU
    # 4 column: indicator whether SKU is allocated
    return sales::Array{Int64,2}
end

## function to sort the SKU pairs after the highest number of coappearances
function SORTPAIRS(Q::Array{Int64,2})
    pairs = Array{Int64,2}(undef,round(Int64,(size(Q,1)*size(Q,1)-size(Q,1))/2),3)
    iter = 1
    for i = 2:size(Q,1)
        for j = 1 : i-1
            pairs[iter,1] = i
            pairs[iter,2] = j
            pairs[iter,3] = Q[i,j]
            iter += 1
        end
    end
    #pairs = sortslices(pairs,dims=1,by=x->x[3],rev=true)
    pairs = pairs[sortperm(pairs[:,3],rev=true),:]
    return pairs::Array{Int64,2}
end

## function to check the number of coappearances in the selected warehouse
## for the selected product
function CURRENTCOAPP(X::Array{Bool,2},
                      warehouse::Int64,
                      sku::Int64,
                      Q::Array{Int64,2})
    sum(X[:,warehouse] .* Q[sku,:])
end

## function to start place the seeds in the greedy seeds heuristic
function GREEDYSEEDSTART!(sales::Array{Int64,2},
                          X::Array{Bool,2},
                          Q::Array{Int64,2},
                          capacity_left::Array{Int64,1})
    for k = 2:size(capacity_left,1)
        # list all SKUs that have not been allocated yet and
        top_sales = sales[sales[:,4] .== 0,:] 
        # discard from the list those SKUs that are not in the top decile
        top_sales = top_sales[1:round(Int64, size(top_sales,1)/10),:]
        # find the SKU from the list with the least coappearances to the allocated SKUs
        coapp = Array{Int64,1}(undef,size(top_sales,1)) .= 0
        for i in 1:size(top_sales,1)
            for j in 1:size(capacity_left,1)-1
                coapp[i] = CURRENTCOAPP(X,j,i,Q)
            end
        end
        # allocate the SKU to the DC
        sales[top_sales[findmin(coapp)[2],2],4] = 1
        X[top_sales[findmin(coapp)[2],1],k] = 1
    end
    return sales::Array{Int64,2},
           X::Array{Bool,2}
end

## function to allocate the SKU with the highest potential
function GREEDYSEEDMAIN!(sales::Array{Int64,2},
                         X::Array{Bool,2},
                         Q::Array{Int64,2},
                         capacity_left::Array{Int64,1})
    for i = 1:size(Q,1)
        if sales[i,4] == 0
            best_allocation = Array{Int64,1}(undef,size(capacity_left,1)) .= 0
            for d = 1:size(capacity_left,1)
                if capacity_left[d] > 1
                    # retrieve the list of SKUs allocated to d so far
                    # calculate the average of coapperances with the SKUs in d
                    best_allocation[d] = CURRENTCOAPP(X,d,i,Q)
                end
            end
            if findmax(best_allocation)[1] > 0
                X[sales[i,1],findmax(best_allocation)[2]] = 1
                capacity_left[findmax(best_allocation)[2]] -= 1
            else
                X[sales[i,1],findmax(capacity_left)[2]] = 1
                capacity_left[findmax(capacity_left)[2]] -= 1
            end
            sales[i,4] = 1
        end
    end
    return sales::Array{Int64,2},
           X::Array{Bool,2},
           capacity_left::Array{Int64,1}
end

## function to allocate the SKUs in the greedy pairs heuristic
function GREEDYPAIRSMAIN!(pairs::Array{Int64,2},
                          X::Array{Bool,2},
                          capacity_left::Array{Int64,1})
    ## Sort the DCs by decreasing capacity
    for i = 1:size(pairs,1)
        ## For each pair of SKUs in sorted list of pairs do
        for j = 1:2
            for d in 1:size(capacity_left,1)
                if capacity_left[d] > 0
                    ## if SKU n is not allocated to DC j then allocate SKU s to DC d
                    if sum(X[pairs[i,j],:]) == 0
                        X[pairs[i,j],d] = 1
                        capacity_left[d] -= 1
                    end
                end
            end
        end
        if sum(X) == size(X,1)
            break
        end
    end
    return X::Array{Bool,2}, 
           capacity_left::Array{Int64,1}
end

## function to fill up the warehouses with the potentially best SKUs
## after each SKU has already been allocated once
function FILLUP(X::Array{Bool,2},
                Q::Array{Int64,2},
                capacity_left::Array{Int64,1})
    for d = 1:size(capacity_left,1)
        while capacity_left[d] > 0
            best_allocation = Array{Int64,1}(undef,size(Q,1)) .= 0
            for i = 1:size(Q,1)
                if X[i,d] == 0
                    a = @view X[:,d]
                    b = @view Q[:,i]
                    best_allocation[i] = dot(a,b)
                end
            end
            if findmax(best_allocation)[1] > 0
                X[findmax(best_allocation)[2],d] = 1
                capacity_left[d] -= 1
            else
                capacity_left[d] = 0
            end
        end
    end
    return X::Array{Bool,2}
end

## function to calculate the number of SKUs that can be allocated
## multiple times
function BESTSELLING_B(skus::Int64, capacity_left::Array{Int64,1})
    floor(Int64,(sum(capacity_left)-skus)/(size(capacity_left,1)-1))
end

## function to initialise the bestselling heuristic
function BESTSELLINGSTART!(sales::Array{Int64,2},
                           capacity_left::Array{Int64,1},
                           B::Int64,
                           X::Array{Bool,2})
    for d = 1:size(capacity_left,1)
        if capacity_left[d] < B
            for i = 1:capacity_left[d]
                X[sales[i,1],d] = 1
                sales[sales[i,2],4] = 1
                capacity_left[d] -= 1
            end
        end
        B = BESTSELLING_B(size(X,1),capacity_left)
    end
    return sales::Array{Int64,2},
           capacity_left::Array{Int64,1},
           B::Int64,
           X::Array{Bool,2}
end

## function to allocate the SKUs that have been sold most in each
## warehouse until the remaining SKUs can only be allocated once
function BESTSELLINGTOP!(sales::Array{Int64,2},
                         capacity_left::Array{Int64,1},
                         B::Int64,
                         X::Array{Bool,2})
    for d = 1:size(capacity_left,1)
        if capacity_left[d] >= B
            for i = 1:B
                X[sales[i,1],d] = 1
                sales[sales[i,2],4] = 1
                capacity_left[d] -= 1
            end
        end
    end
    return sales::Array{Int64,2},
           capacity_left::Array{Int64,1},
           X::Array{Bool,2}
end


