# Function to calculate the Coappearance Matrix Q
function COAPPEARENCE(
    trans::SparseMatrixCSC{Bool, Int64},
    sku_weight::Vector{<:Real}
    )
    Q = trans' * trans
    return Q
end

# function to set all entries on the principle diagonal to zero
function CLEANPRINCIPLE!(
    Q::AbstractMatrix{<:Real}
    )
    for i in axes(Q,1)
        Q[i,i] = 0
    end
    if Q isa SparseMatrixCSC
        dropzeros!(Q)
    end
end

# Functions to evaluate the number of parcels necessary to fullfil all orders
## COMBINEWAREHOUSES: create all possible warehouse combinations 
function COMBINEWAREHOUSES(
    capacity::Array{Int64,1}
    )
    combination = [[k] for k in axes(capacity,1)]
    combinations_X = combinations(combination)
    combination = [c for c in combinations_X]
    return combination::Array{Array{Array{Int64,1},1},1}
end

## PARCELSSEND_WEIGHT: shows the number of parcels necessary plus the weight for each warehouse dispatch
## Attention: does only work with 2 Warehouses!!
function PARCELSSEND_WEIGHT(
    trans::SparseMatrixCSC{Bool, Int64}, 
    X::Array{Bool,2}, 
    capacity::Array{Int64,1}, 
    combination::Array{Array{Array{Int64,1},1},1},
    max_warehouse_1::Bool
    )
    if sum(capacity) == size(trans,2)
        parcel = trans * X
        parcel = (parcel .> 0)
        parcel = sum(parcel, dims=2)
    else
        warehouse_combination = zeros(Int64, size(X,1), size(combination,1))
        parcels_out = zeros(size(X,2), size(combination,1))
        for i in axes(combination,1)
            for j in axes(combination[i],1)
                for k in axes(X,1)
                    if X[k,combination[i][j]] == [1]
                        warehouse_combination[k,i] = 1
                    end
                end
                parcels_out[combination[i][j][1],i] = 1
            end
        end
        warehouse_combination = dropzeros(sparse(warehouse_combination))
        parcel_evalu = trans * warehouse_combination
        parcel_check = sum(trans, dims = 2)
        parcel = zeros(Int64,size(X,2))
        for i in axes(parcel_evalu,1)
            if max_warehouse_1 == true
                for j in [2,1,3] #axes(parcel_evalu,2)
                    if parcel_evalu[i,j] == parcel_check[i]
                        parcel[:] .+= parcels_out[:,j]
                        break
                    else
                        if j == maximum([2,1,3])
                            parcel[:] .+= parcels_out[:,j]
                        end
                    end
                end
            else
                for j in axes(parcel_evalu,2)
                    if parcel_evalu[i,j] == parcel_check[i]
                        parcel[:] .+= parcels_out[:,j]
                        break
                    else
                        if j == maximum(axes(parcel_evalu,2))
                            parcel[:] .+= parcels_out[:,j]
                        end
                    end
                end
            end
        end
        parcel = transpose(parcel)
    end
    split = sum(parcel) - size(trans,1)
    return split, parcel
end

## PARCELSSEND: number of parcels due to ordersplitting necessary to fulfill all orders
function PARCELSSEND(
    trans::SparseMatrixCSC{Bool, Int64}, 
    X::Array{Bool,2}, 
    capacity::Array{Int64,1}, 
    combination::Array{Array{Array{Int64,1},1},1}
    )
    if !any(iszero, sum(X, dims=2))
        if sum(capacity) == size(trans,2)
            parcel = trans * X
            parcel = sum(parcel .> 0)
        else
            warehouse_combination = zeros(Int64, size(X,1), size(combination,1))
            parcels_out = zeros(size(combination,1))
            for i in axes(combination,1)
                for j in axes(combination[i],1)
                    for k in axes(X,1)
                        if X[k,combination[i][j]] == [1]
                            warehouse_combination[k,i] = 1
                        end
                    end
                    parcels_out[i] += 1
                end
            end
            warehouse_combination = dropzeros(sparse(warehouse_combination))
            parcel_evalu = trans * warehouse_combination
            parcel_check = sum(trans, dims = 2)
            parcel = 0
            for i in axes(parcel_evalu,1)
                for j in axes(parcel_evalu,2)
                    if parcel_evalu[i,j] == parcel_check[i]
                        parcel += parcels_out[j]
                        break
                    end
                    if j == size(parcel_evalu,2)
                        if parcel_evalu[i,j] != parcel_check[i]
                            error("Order was not fully dispatched!")
                        end
                    end
                end
            end
        end
        parcel -= size(trans,1)
    else
        show(X)
        error("Not all SKUs are allocated!")
    end
    return parcel
end

# Functions for the random allocation of SKUs to warehouses
## RANDOMALLOCONCE: allocate SKUs randomly in case each SKU can only be assigned once
function RANDOMALLOCONCE(
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real},
    )
    cap_left = copy(capacity)
    X = zeros(Bool, length(sku_weight), length(capacity))
    
    # Pre-calculate warehouse indices
    warehouse_indices = collect(1:length(capacity))
    
    for j in randperm(length(sku_weight))
        while sum(view(X, j, :)) == 0  # Use view instead of full slice
            randomnumber = rand(warehouse_indices)
            if cap_left[randomnumber] > 0
                X[j,randomnumber] = true
                cap_left[randomnumber] -= sku_weight[j]
            end
        end
    end
    
    any(iszero, sum(X, dims=2)) && error("Not all SKUs are allocated.")
    return X
end

## RANDOMALLOCMULTI: allocate SKUs randomly in case each SKU can be allocated multiple times
function RANDOMALLOCMULTI(
    trans::SparseMatrixCSC{Bool, Int64},
    capacity::Vector{<:Real},
    sku_weight::Vector{<:Real}
    )
    X = RANDOMALLOCONCE(capacity,sku_weight::Vector{<:Real})
    cap_used = sum(X.*sku_weight,dims=1)
    for d in axes(capacity,1)
        while cap_used[d] < capacity[d]
            randomproduct = rand(1:size(trans,2))
            if X[randomproduct,d] == 0
                X[randomproduct,d] = 1
                cap_used[d] += sku_weight[randomproduct]
            end
        end
    end
    return X::Array{Bool,2}
end

## RANDOMBENCH: benchmark function to evaluate multiple random allocations
function RANDOMBENCH(
    trans::SparseMatrixCSC{Bool, Int64}, 
    capacity::Array{Int64,1}, 
    iterations::Int64,
    sku_weight::Vector{<:Real},
    combination::Array{Array{Array{Int64,1},1},1}
    )
    randomcollection = zeros(Int64, iterations)
    Threads.@threads for r = 1:iterations
        W = RANDOMALLOCMULTI(trans, capacity, sku_weight)
        randomcollection[r] = PARCELSSEND(trans, W, capacity, combination)
    end
    return round(mean(randomcollection))  # More efficient than sum/length
end

## CHECKCAPACITY: function to test whether the input capacity is viable for the transactional data
function CHECKCAPACITY(
    capacity::Array{Int64,1},
    sku_weight::Vector{<:Real},
    )
    valid = true
    if sum(capacity) < sum(sku_weight)
        error("Abort due to insufficient capacity!")
        valid = false
    end
    return valid::Bool
end