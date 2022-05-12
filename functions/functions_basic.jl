# Function to calculate the Coappearance Matrix Q
function COAPPEARENCE(trans::Matrix{Int64})
    t_sparse = dropzeros(sparse(trans))
    Q = t_sparse'*t_sparse
    for i = 1:size(Q,1)
        Q[i,i] = 0
    end
    Q = Matrix(Q)
    Q = Int.(Q)
    return Q
end

function COAPPEARENCE(trans::SparseMatrixCSC{Float64, Int64})
    Q = trans'*trans
    for i = 1:size(Q,1)
        Q[i,i] = 0
    end
    Q = Matrix(Q)
    Q = Int.(Q)
    return Q
end

# Functions to evaluate the number of parcels necessary to fullfil all orders
## COMBINEWAREHOUSES: create all possible warehouse combinations 
function COMBINEWAREHOUSES(capacity::Array{Int64,1})
    combination = [[k] for k = 1:size(capacity,1)]
    combinations_X = combinations(combination)
    combination = [c for c in combinations_X]
    return combination::Array{Array{Array{Int64,1},1},1}
end

## PARCELSSEND: number of parcels necessary to fulfill all orders
function PARCELSSEND(trans::SparseMatrixCSC{Float64, Int64}, 
                     X::Array{Int64,2}, 
                     capacity::Array{Int64,1}, 
                     combination::Array{Array{Array{Int64,1},1},1})
    trans = trans
    if sum(capacity) == size(trans,2)
        parcel = trans * X
        for j = 1:size(parcel,1)
            for k = 1:size(parcel,2)
                if parcel[j,k] > 0
                    parcel[j,k] = 1
                else
                    parcel[j,k] = 0
                end
            end
        end
        parcel = round(Int64, sum(parcel))
    else
        warehouse_combination = Array{Int64,2}(undef,size(X,1),size(combination,1)) .= 0
        parcels_out = zeros(size(combination,1))
        for i = 1:size(combination,1)
            for j = 1:size(combination[i],1)
                for k = 1:size(X,1)
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
        for i in 1:size(parcel_evalu,1)
            for j in 1:size(parcel_evalu,2)
                if parcel_evalu[i,j] == parcel_check[i]
                    parcel += parcels_out[j]
                    break
                end
            end
        end
    end
    return parcel
end

# Functions for the random allocation of SKUs to warehouses
## RANDOMALLOCONCE: allocate SKUs randomly in case each SKU can only be assigned once
function RANDOMALLOCONCE(trans::SparseMatrixCSC{Float64, Int64},
                         capacity::Array{Int64,1})
    X = Array{Bool,2}(undef,size(trans,2),size(capacity,1))
    X .= 0
    for j in 1:size(X,1)
        while sum(X[j,:]) == 0
            randomnumber = rand(1:size(capacity,1))
            if sum(X[:,randomnumber]) < capacity[randomnumber]
                X[j,randomnumber] = 1
            end
        end
    end
    X = X[shuffle(1:size(X,1)),:]
    return X::Array{Bool,2}
end

## RANDOMALLOCMULTI: allocate SKUs randomly in case each SKU can be allocated multiple times
function RANDOMALLOCMULTI(trans::SparseMatrixCSC{Float64, Int64},
                          capacity::Array{Int64,1})
    X = RANDOMALLOCONCE(trans,capacity)
    for d = 1:size(capacity,1)
        while sum(X[:,d]) < sum(capacity[d])
            randomproduct = rand(1:size(trans,2))
            if X[randomproduct,d] == 0
                X[randomproduct,d] = 1
            end
        end
    end
    return X::Array{Bool,2}
end

## RANDOMBENCH: benchmark function to evaluate multiple random allocations
function RANDOMBENCH(trans::SparseMatrixCSC{Float64, Int64}, 
                     capacity::Array{Int64,1}, 
                     iterations::Int64, 
                     combination::Array{Array{Array{Int64,1},1},1})
    randomcollection = Array{Int64,1}(undef,iterations) .= 0
    Threads.@threads for r = 1:iterations
        W = RANDOMALLOCMULTI(trans,capacity)
        W = convert(Matrix{Int64},W)
        parcel = PARCELSSEND(trans,W,capacity,combination)
        Threads.lock(ren_lock) do
            randomcollection[r] = parcel
        end
    end
    average_parcels = round(sum(randomcollection)/iterations,digits=0)
    return average_parcels::Float64
end

## CHECKCAPACITY: function to test whether the input capacity is viable for the transactional data
function CHECKCAPACITY(trans,
                       capacity::Array{Int64,1})
    valid = 1
    if sum(capacity) < size(trans,2)
        error("Abort due to insufficient capacity!")
        valid = 0
    end
    return valid::Int64
end

# function to apply a pair-wise exchange local search on the allocation
# of all other heuristics
function LOCALSEARCH(W::Array{Int64,2},
                     Q::Array{Int64,2})
    # lock to avoid racing conditions during the local search
    ren_lock = ReentrantLock()
    iteration = 1
    coapp_sort = vec(sum(Q,dims = 2))
    coapp_sort = sortperm(coapp_sort,rev=true)
    while iteration > 0
        iteration = 0
        for i in coapp_sort
            for j in coapp_sort
                if i != j
                    for k = 2:size(W,2)
                        for g = 1:size(W,2)-1
                            if W[i,k] == 1 && W[j,g] == 1 && W[i,g] == 0 && W[j,k] == 0
                                iteration += POTENTIAL!(i,j,k,g,W,Q)
                            end
                        end
                    end
                end
            end
        end
    end
   return W::Array{Int64,2}
end

# function to calculate the potential improvement during
# the local search procedure
function POTENTIAL!(i::Int64,
                        j::Int64,
                        k::Int64,
                        g::Int64,
                        X::Array{Int64,2},
                        Q::Array{Int64,2})
    potential = 0
    impr = 0
    potential += dot(@view(X[:,g]),@view(Q[i,:]))
    potential -= dot(@view(X[:,k]),@view(Q[i,:]))
    potential += dot(@view(X[:,k]),@view(Q[j,:]))
    potential -= dot(@view(X[:,g]),@view(Q[j,:]))
    if potential > 0
        impr = 1
        X[i,k] = 0
        X[j,g] = 0
        X[i,g] = 1
        X[j,k] = 1
    end
    return impr::Int64
end