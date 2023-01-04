## function to sort the SKUs after the number of sales
function SORTSALES(trans::SparseMatrixCSC{Bool, Int64}, sku_weight::Vector{<:Real})
    sales = Array{Int64,2}(undef,size(trans,2),4) .= 0
    sales[:,3] = round.(Int64, sum(trans,dims=1)./transpose(sku_weight))
    sales[:,4] .= 0
    for i in axes(trans,2)
        sales[i,1] = i
    end
    sales = sortslices(sales,dims=1,by=x->x[3],rev=true)
    for i in axes(trans,2)
        sales[i,2] = i
    end
    # 1 column: index of SKU in allocation matrix
    # 2 column: index of SKU in sorted sales
    # 3 column: number of total sold SKU
    # 4 column: indicator whether SKU is allocated
    return sales::Matrix{Int64}
end

## function to sort the SKU pairs after the highest number of coappearances
function SORTPAIRS(Q::Matrix{<:Real}, sku_weight::Vector{<:Real})
    pairs = zeros(Int64,round(Int64,(size(Q,1)*size(Q,1)-size(Q,1))/2),3)
    iter = 1
    for i = 2:size(Q,1)
        for j = 1 : i-1
            pairs[iter,1] = i
            pairs[iter,2] = j
            pairs[iter,3] = round(Int64, Q[i,j]./(sku_weight[i]+sku_weight[j])/2)
            iter += 1
        end
    end
    pairs = pairs[sortperm(@view(pairs[:,3]),rev=true),:]
    return pairs::Array{Int64,2}
end

## function to check the number of coappearances in the selected warehouse
## for the selected product
function CURRENTCOAPP(X::Matrix{Bool},
                      warehouse::Int64,
                      sku::Int64,
                      Q::Matrix{<:Real})
    sum(X[:,warehouse] .* Q[sku,:])
end

## function to start place the seeds in the greedy seeds heuristic
function GREEDYSEEDSTART!(sales::Matrix{Int64},
                          X::Matrix{Bool},
                          Q::Matrix{<:Real},
                          capacity_left::Vector{<:Real},
                          sku_weight::Vector{<:Real})
    for k = 2:size(capacity_left,1)
        if capacity_left[k] > 0
            # list all SKUs that have not been allocated yet and
            top_sales = sales[sales[:,4] .== 0,:] 
            # discard from the list those SKUs that are not in the top decile
            top_sales = top_sales[1:round(Int64, size(top_sales,1)/10),:]
            # find the SKU from the list with the least coappearances to the allocated SKUs
            coapp = zeros(Float64,size(top_sales,1))
            for i in 1:size(top_sales,1)
                for j in 1:size(capacity_left,1)-1
                    coapp[i] = CURRENTCOAPP(X,j,i,Q)/sku_weight[i]
                end
            end
            # allocate the SKU to the DC
            sales[top_sales[findmin(coapp)[2],2],4] = 1
            X[top_sales[findmin(coapp)[2],1],k] = 1
            capacity_left[k] -= sku_weight[top_sales[findmin(coapp)[2],1]]
        end
    end
end

## function to allocate the SKU with the highest potential
function GREEDYSEEDMAIN!(sales::Matrix{Int64},
                         X::Matrix{Bool},
                         Q::Matrix{<:Real},
                         capacity_left::Vector{<:Real},
                         sku_weight::Vector{<:Real})
    for i in axes(Q,1)
        if sales[i,4] == 0
            best_allocation = zeros(Float64,size(capacity_left,1))
            for d in axes(capacity_left,1)
                if capacity_left[d] >= sku_weight[i]
                    # retrieve the list of SKUs allocated to d so far
                    # calculate the average of coapperances with the SKUs in d
                    best_allocation[d] = CURRENTCOAPP(X,d,i,Q)/sku_weight[i]
                end
            end
            if findmax(best_allocation)[1] > 0
                X[sales[i,1],findmax(best_allocation)[2]] = 1
                capacity_left[findmax(best_allocation)[2]] -= sku_weight[sales[i,1]]
            else
                X[sales[i,1],findmax(capacity_left)[2]] = 1
                capacity_left[findmax(capacity_left)[2]] -= sku_weight[sales[i,1]]
            end
            sales[i,4] = 1
        end
    end
end


function GREEDYPAIRSMAIN!(pairs::Matrix{Int64},
    X::Matrix{Bool},
    capacity_left::Vector{<:Real},
    sku_weight::Vector{<:Real})
    for i in axes(pairs,1)
        ## For each pair of SKUs in sorted list of pairs do
        for j = 1:2
            for d in 1:size(capacity_left,1)
                if capacity_left[d] >= sku_weight[pairs[i,j]]
                    ## if SKU n is not allocated to DC j then allocate SKU s to DC d
                    if sum(X[pairs[i,j],:]) == 0
                        X[pairs[i,j],d] = 1
                        capacity_left[d] -= sku_weight[pairs[i,j]]
                        break
                    end
                end
            end
        end
        if sum(X) == size(X,1)
            break
        end
    end
end

# function to allocate the SKUs in the main part of the Greedy Orders Heuristic
function GREEDYORDERSMAIN!(
    order_sort::Vector{Int64},
    trans::SparseMatrixCSC{Bool,Int64},
    capacity_left::Vector{<:Real},
    X::Matrix{Bool},
    sku_weight::Vector{<:Real})
    for order in order_sort
        for sku in axes(trans,2)
            # Determine the warehouse d with the lowest index that has one unallocated space
            if trans[order,sku] == 1
                for d = 1:length(capacity_left)
                    if capacity_left[d] >= sku_weight[sku]
                        ## if the SKU is not allocated to the warehouses allocate it
                        if sum(X[sku,:]) == 0
                            X[sku,d] = 1
                            capacity_left[d] -= sku_weight[sku]
                        end
                    end
                end
            end
        end
    end
    skus_assigned = sum(X,dims = 2)
    for sku in axes(trans,2)
        if skus_assigned[sku] == 0
            while sum(X[sku,:]) == 0
                k = rand(1:length(capacity_left))
                if capacity_left[k] >= sku_weight[sku]
                    X[sku,k] = 1
                    capacity_left[k] -= sku_weight[sku]
                end
            end
        end
    end
end

## function to calculate the number of SKUs that can be allocated
## multiple times
function BESTSELLING_B(
    capacity_left::Vector{<:Real},
    sku_weights::Vector{<:Real},
)
    floor(Int64,(sum(capacity_left)-sum(sku_weights))/(size(capacity_left,1)-1))
end

## function to initialise the bestselling heuristic
function BESTSELLINGSTART!(
    sales::Matrix{Int64},
    capacity_left::Vector{<:Real},
    B::Int64,
    X::Matrix{Bool},
    sku_weight::Vector{<:Real},
)
    for d in axes(capacity_left,1)
        if capacity_left[d] < B
            i = 1
            while capacity_left[d] >= sku_weight[sales[i,1]]
                X[sales[i,1],d] = 1
                sales[sales[i,2],4] = 1
                capacity_left[d] -= sku_weight[sales[i,1]]
                i += 1
            end
        end
        B = BESTSELLING_B(capacity_left,sku_weight)
    end
end

## function to allocate the SKUs that have been sold most in each
## warehouse until the remaining SKUs can only be allocated once
function BESTSELLINGTOP!(sales::Matrix{Int64},
                         capacity_left::Vector{<:Real},
                         capacity::Vector{<:Real},
                         B::Int64,
                         X::Matrix{Bool}, 
                         sku_weight::Vector{<:Real})
    for d in axes(capacity_left,1)
        i = 1
        while capacity_left[d] > capacity[d] - B && i <= size(sales,1)
            X[sales[i,1],d] = 1
            sales[sales[i,2],4] = 1
            capacity_left[d] -= sku_weight[sales[i,1]]
            i += 1
        end
    end
end


