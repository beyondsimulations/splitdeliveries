# RANDOMTRANS: function that generates random transactions with dependencies between
#              SKUs based on the values of max_groupsize and max_dependence
function RANDOMTRANS(skus::Int64,
                     orders::Int64,
                     max_group_size::Int64,
                     min_strength::Float64,
                     max_strength::Float64,
                     group_link::Float64,
                     ind_chance::Float64)
    C = Matrix{Float64}(undef,skus,skus) .= 0.0
    D = Matrix{Int64}(undef,skus,skus)   .= 0
    if max_strength > 0.0
        groups = 0
        while groups < skus
            test_correlation = sum(D,dims=1)
            free_skus = findall(==(0),test_correlation)
            group_size = rand(1:max_group_size)
            if group_size <= length(free_skus)
                members = Vector{Int64}(undef,group_size)
                for i = 1:group_size
                    members[i] = rand(1:length(free_skus))
                end
                for i in randperm(skus)[1:rand(1:ceil(Int64,group_size*group_link))]
                    j = rand(1:skus)
                    strength = rand(Uniform(min_strength,max_strength),1)[1]
                    C[i,j] = strength
                    C[j,i] = strength
                    D[i,j] = 1
                    D[j,i] = 1
                end
                for i = 1:group_size
                    for j = 1:group_size
                        if i != j
                            strength = rand(Uniform(min_strength,max_strength),1)[1]
                            local_i = members[i]
                            local_j = members[j]
                            D[local_i,local_j] = 1
                            D[local_j,local_i] = 1
                            C[local_i,local_j] = strength
                            C[local_j,local_i] = strength
                        end
                    end
                end
                groups += group_size
            else
                break
                groups = max_groupmembers
            end
        end
    end
    #transactions = Matrix{Int64}(undef,orders,skus) .= 0
    transactions = spzeros(orders,skus)
    for i = 1:orders
        already_allocated = 0
        skus_order = 1 + floor(abs(rand(Normal(0,3))))
        while already_allocated < skus_order 
            new_sku = skus + 1
            while new_sku > skus
                new_sku = 1+ floor(Int64, abs(rand(Normal(0,skus/2))))
            end
            already_allocated += 1
            transactions[i,new_sku] = 1
            if rand() > ind_chance
                for j in randperm(skus)
                    if D[new_sku,j] == 1
                        if transactions[i,j] == 0
                            if already_allocated < skus_order
                                if rand() < C[new_sku,j]
                                    transactions[i,j] = 1
                                    already_allocated += 1
                                end
                            else
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    return transactions
end
