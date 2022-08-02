# function used to perform the chi square test of independence upon 
# the coappearance matrix Q
function HYOPTHESISTEST!(dep::Matrix{<:Real},
                        Q::Matrix{<:Real},
                        I::Int64,
                        J::Int64,
                        sig::Float64,
                        sum_cond_sku::Vector{<:Real})
    ## Initialise chi square test
    M = convert(Int64,((I^2-I)/2))
    ## Determine the acceptance level of the test
    accept = cquantile(Chisq(1), sig/M)
    ## Fill the arrays with the results of the chi-square test
    DEPTEST!(dep,Q,sum_cond_sku,J,accept)
end

function DEPTEST!(dep::Matrix{<:Real},
                  Q::Matrix{<:Real},
                  sum_cond_sku::Vector{<:Real},
                  J::Int64,
                  accept::Float64)
    @inbounds @simd for j = 2:size(Q,1)
        @inbounds for i = 1:j-1
            chi_values!(dep,Q,sum_cond_sku,J,accept,i,j)
        end
    end
end

function chi_part(a::Real,
                  b::Real)
    (a-b)^2/b
end

function chi_values!(dep::Matrix{<:Real},
                     Q::Matrix{<:Real}, 
                     sum_cond_sku::Vector{<:Real}, 
                     J::Int64,
                     accept::Float64,
                     i::Int64,
                     j::Int64)
    chi = 0.0
    independent = INDEPENDENT(sum_cond_sku,J,i,j)
    if Q[i,j] > independent
        chi_nr = J - sum_cond_sku[i]
        chi_nd = J - sum_cond_sku[j]
        chi_yn = sum_cond_sku[i] - Q[i,j]
        chi_ny = sum_cond_sku[j] - Q[i,j]
        chi_nn = chi_nd - chi_yn
        ind_yn = (sum_cond_sku[i] * chi_nd)/J
        ind_ny = (sum_cond_sku[j] * chi_nr)/J
        ind_nn = (chi_nr * chi_nd)/J
        chi += chi_part(Q[i,j],independent)
        chi += chi_part(chi_yn,ind_yn)
        chi += chi_part(chi_ny,ind_ny)
        chi += chi_part(chi_nn,ind_nn)
        if chi > accept
            dep[i,j] = dep[j,i] = Q[i,j] - independent
        end
    end
end

function fish_values!(dep::Matrix{<:Real},
                     Q::Matrix{<:Real}, 
                     sum_cond_sku::Vector{<:Real}, 
                     J::Int64,
                     accept::Float64,
                     i::Int64,
                     j::Int64)
    chi = 0.0
    independent = INDEPENDENT(sum_cond_sku,J,i,j)
    if Q[i,j] > independent && independent > 10
        chi_nr = J - sum_cond_sku[i]
        chi_nd = J - sum_cond_sku[j]
        chi_yn = sum_cond_sku[i] - Q[i,j]
        chi_ny = sum_cond_sku[j] - Q[i,j]
        chi_nn = chi_nd - chi_yn
        ind_yn = (sum_cond_sku[i] * chi_nd)/J
        ind_ny = (sum_cond_sku[j] * chi_nr)/J
        ind_nn = (chi_nr * chi_nd)/J
        chi += chi_part(Q[i,j],independent)
        chi += chi_part(chi_yn,ind_yn)
        chi += chi_part(chi_ny,ind_ny)
        chi += chi_part(chi_nn,ind_nn)
        if chi > accept
            dep[i,j] = dep[j,i] = Q[i,j] - independent
        end
    end
end

function INDEPENDENT(sum_cond_sku::Vector{<:Real},
                     J::Int64,
                     i::Int64,
                     j::Int64)
    (sum_cond_sku[i] * sum_cond_sku[j])/J
end

# function to calculate the weight of each warehouse. It shows us the density of 
# the independent coappearances in each warehouse if we were to allocate
# all SKUs according to the highest independent coappearances.
function WHWEIGHT(capacity::Vector{Int64},
                  sum_nor::Vector{<:Real})
    sort_sum_nor = sort!(copy(sum_nor), rev = true)
    sort_sum_nor = vcat(sort_sum_nor,zeros(Float64,sum(capacity)-I))
    weight = zeros(Float64,size(capacity))
    weight[1] = sum(sort_sum_nor[1:capacity[1]])/sum(sort_sum_nor)
    for k = 2:size(capacity,1)
        weight[k] = sum(sort_sum_nor[sum(capacity[1:k-1]):sum(capacity[1:k])])/
                    sum(sort_sum_nor)
    end
    return weight::Vector{Float64}
end

# function to find the largest warehouse that still has
# space left
function WHSPACE(cap_left::Vector{Int64},
                 space::Int64)
    k_ind = 0
    for k in 1:size(cap_left,1)
        if cap_left[k] >= space
            k_ind = k
            break
        end
    end
    return k_ind::Int64
end

# function to select the best SKU for an allocation
function SELECTIK(sum_dep::Vector{<:Real},
                  sum_nor::Vector{<:Real},
                  weight::Vector{Float64},
                  cap_left::Vector{Int64},
                  X::Array{Bool,2},
                  dep::Matrix{<:Real},
                  allocated::Vector{Bool})
    i       = findmax(sum_nor)[2]
    if sum(@view(X[i,:])) > 0
        for j = 1:size(X,1)
            if allocated[j] == 0
                i = j
                break
            end
        end
    end
    k_ind   = WHSPACE(cap_left,1)
    k_max   = findmax(cap_left)[2]
    pot_dep = WHPOTDEP(cap_left,i,X,dep)
    k_dep   = findmax(pot_dep)[2]
    if pot_dep[k_dep] > 0 && k_ind != k_ind
        sum_dep[i] + sum_nor[i] * weight[k_dep] > 
        sum_nor[i] * weight[k_ind]
            k = k_dep
        elseif k_max != k_ind && sum_dep[i] + sum_nor[i] * weight[k_max] >
            sum_nor[i] * weight[k_ind]
                k = k_max
        else
            k = k_ind
    end
    return i::Int64,
           k::Int64
end

# function to check for each warehouse with free space whether 
# SKU i has significant dependencies to other already allocated SKUs
function WHPOTDEP(cap_left::Vector{Int64},
                  i::Int64,
                  X::Array{Bool,2},
                  dep::Matrix{<:Real})
    pot_dep = zeros(Float64,size(cap_left,1))
    for k in 1:size(cap_left,1)
        if cap_left[k] > 0
            pot_dep[k] = CALCVAL(X,dep,i,k)
        end
    end
    return pot_dep::Vector{Float64}
end

# function to remove every so far assigned dependent SKU-pair from 
# the coappearance matrix dep to prevent the allocation bias described 
# in our article. 
function REMOVEALLOC!(X::Array{Bool,2},
                      Q::Matrix{<:Real},
                      dep::Matrix{<:Real})
    for b = 2:size(Q,1)
        for a = 1:b
            for k = 1:size(X,2)
                if  X[a,k] == 1 && X[b,k] == 1 && dep[a,b] > 0
                    dep[a,b] = dep[b,a] = Q[a,b] - dep[a,b]
                else
                    dep[a,b] = dep[b,a] = Q[a,b]
                end
            end
        end
    end
end

# function to check for all unallocated SKUs whether they have positive 
# dependencies to the SKUs in the warehouse k the last SKU was allocated 
# to. If so, check whether the dependencies are expected to dominate the 
# independent coapperances. If yes, allocate the corresponding SKUs to 
# the warehouse k.
function ADDDEPENDENT!(X::Matrix{Bool},
                       Q::Matrix{<:Real},
                       cap_left::Vector{Int64},
                       k::Int64,
                       dep::Matrix{<:Real},
                       sum_dep::Vector{<:Real},
                       sum_nor::Vector{<:Real},
                       state_dep::Matrix{<:Real},
                       state_nor::Matrix{<:Real},
                       allocated::Vector{Bool})
    add = 1
    pot_dep = zeros(Float64,size(X,1))
    pot_nor = zeros(Float64,size(X,1))
    while add == 1 && cap_left[k] > 0 && sum(X) < size(X,1)
        add = 0
        FINDDEP!(X,k,state_dep,state_nor,pot_dep,pot_nor,allocated)
        i = findmax(pot_dep)[2]
        if findmax(pot_dep)[1] > 0
            if pot_dep[i] >= findmax(pot_nor)[1]
                ALLOCATEONE!(X,dep,Q,sum_dep,sum_nor,state_dep,
                             state_nor,cap_left,allocated,i,k)
                add = 1
            end
        end
    end
end

# function to check the dependencies to already allocated SKUs
function FINDDEP!(X::Array{Bool,2},
                  k::Int64,
                  state_dep::Matrix{<:Real},
                  state_nor::Matrix{<:Real},
                  pot_dep::Vector{<:Real},
                  pot_nor::Vector{<:Real},
                  allocated::Vector{Bool})
    @fastmath @simd for j in 1:size(X,1)
        if allocated[j] == 0
            pot_dep[j]  = state_dep[j,k] 
            pot_nor[j]  = state_nor[j,k] 
            if pot_dep[j] > 0
                pot_dep[j] += pot_nor[j]
            end
        else
            pot_dep[j]  = 0
            pot_nor[j]  = 0
        end
    end
end

# function to allocate a selcted product to a selected warehouse
function ALLOCATEONE!(X::Array{Bool,2},
                      dep::Matrix{<:Real},
                      Q::Matrix{<:Real},
                      sum_dep::Vector{<:Real},
                      sum_nor::Vector{<:Real},
                      state_dep::Matrix{<:Real},
                      state_nor::Matrix{<:Real},
                      cap_left::Vector{Int64},
                      allocated::Vector{Bool},
                      i::Int64,
                      k::Int64)
    if sum(X[i,:]) == 0
        X[i,k] = 1
        sum_dep[i] = 0
        sum_nor[i] = 0
        cap_left[k] -= 1
        allocated[i] = true
        state_dep[i,k] = 0
        state_nor[i,k] = 0
        @fastmath for j in 1:size(X,1)
            if allocated[j] == 0
                state_dep[j,k]  += dep[j,i]
                state_nor[j,k]  += Q[j,i] - dep[j,i]
            end
        end
    else
        error("This product is already allocated!")
    end
end

## check whether all warehouse except the last one are already full. If
## that is the case just allocate the remaining SKUs yet not allocated there.
function FILLLAST!(X::Array{Bool,2},
                   cap_left::Vector{Int64},
                   allocated::Vector{Bool})
    if sum(cap_left[1:size(cap_left,1)-1]) == 0
        for i = 1:length(allocated)
            if allocated[i] == false
                X[i,size(cap_left,1)] = 1
                cap_left[size(cap_left,1)] -= 1
            end
        end
    end
end


function CALCVAL(X::Matrix{Bool},
                 T::Matrix{<:Real},
                 i::Int64,
                 k::Int64)
    out = 0
    @avxt for y = 1:size(X,1)
        out += X[y,k] * T[y,i]
    end
    return out
end

# function to allocate the SKUs with the highest potential allocation 
# value to each warehouse with leftover storage space until it is full
function FILLUP!(X::Array{Bool,2},
                 Q::Matrix{<:Real},
                 capacity_left::Vector{Int64})
    state = Matrix{Float64}(undef,size(X,1),size(X,2)) .= 0
    for i = 1:size(X,1)
        for d = 1:size(X,2)
            state[i,d] = CALCVAL(X,Q,i,d)
        end
    end
    for d = 1:size(capacity_left,1)
        while capacity_left[d] > 0
            best_allocation = Array{Float64,1}(undef,size(Q,1)) .= 0
            for i = 1:size(Q,1)
                if X[i,d] == 0
                    best_allocation[i] = state[i,d]
                end
            end
            best = findmax(best_allocation)
            if best[1] > 0
                X[best[2],d] = 1
                for j in 1:size(X,1)
                    state[j,d] += Q[j,best[2]]
                end
                capacity_left[d] -= 1
            else
                capacity_left[d] = 0
            end
        end
    end
end

# function to apply a pair-wise exchange local search on the allocation
# of the CHI heuristic
function LOCALSEARCHCHI!(trans::SparseMatrixCSC{Bool,Int64},
                         X::Matrix{Bool}, 
                         Q::Matrix{<:Real},
                         capacity::Vector{Int64},
                         log_results::Bool,
                         ls::Int64,
                         max_ls::Int64)
    coapp_sort = vec(sum(Q,dims = 1))
    coapp_sort = sortperm(coapp_sort,rev=true)
    combination = COMBINEWAREHOUSES(capacity)
    state = zeros(Int32,size(X,1),size(X,2))
    CURRENTSTATE!(X,Q,state)
    impro_bef = PARCELSSEND(trans,X,capacity,combination)
    X_backup = zeros(Bool,size(X,1),size(X,2))
    log_results == true ? print("\n  Iter: 0 - parcels: ",impro_bef,"\n") : nothing
    while ls < max_ls
        ls += 1
        X_backup .= X
        SEARCHLOOP!(X,Q,coapp_sort,state)
        impro_now = PARCELSSEND(trans,X,capacity,combination)
        log_results == true ? print("  Iter: ",ls," - parcels: ",impro_now,"\n") : nothing
        if impro_now < impro_bef
            impro_bef = impro_now
        else
            X .= X_backup
            break
        end
    end
    return ls
end

function POTENTIAL(state::Matrix{<:Real},
                   i::Int64,
                   j::Int64,
                   k::Int64,
                   g::Int64)
    state[i,g] - state[i,k] + state[j,k] - state[j,g]
end

function CURRENTSTATE!(X::Matrix{Bool},
                       Q::Matrix{<:Real},
                       state::Matrix{<:Real})
    @inbounds for k = 1:size(X,2)
        @inbounds for i = 1:size(X,1)
            state[i,k] += CALCVAL(X,Q,i,k)
        end
    end
end

function REFRESHSTATE!(state::Matrix{<:Real},
                       Q::Matrix{<:Real},
                       k::Int64,
                       g::Int64,
                       i::Int64,
                       j::Int64)
    @avxt for y in 1:size(Q,1)
        state[y,k]  += Q[y,j] - Q[y,i]
        state[y,g]  += Q[y,i] - Q[y,j]
    end
end

function SEARCHLOOP!(X::Matrix{Bool},
                     Q::Matrix{<:Real},
                     coapp_sort::Vector{Int64},
                     state::Matrix{<:Real})
    @fastmath begin
        @inbounds for g = 2:size(X,2)
            @inbounds for k = 1:size(X,2)-1
                @inbounds for i in coapp_sort
                    if X[i,k] == 1 && X[i,g] == 0
                        @inbounds for j in coapp_sort
                            if X[j,g] == 1 && X[j,k] == 0
                                if POTENTIAL(state,i,j,k,g) > 0
                                    X[i,k] = X[j,g] = false
                                    X[i,g] = X[j,k] = true
                                    REFRESHSTATE!(state,Q,k,g,i,j)
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end    