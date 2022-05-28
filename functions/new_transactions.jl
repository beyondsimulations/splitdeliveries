# RANDOMTRANS: function that generates random transactions with dependencies between
#              SKUs based on the values of max_groupsize and max_dependence
function RANDOMTRANS(skus::Int64,
                     orders::Int64,
                     max_group_size::Int64,
                     min_dependence::Float64,
                     max_dependence::Float64,
                     group_link::Float64,
                     ind_chance::Float64,
                     one_direction::Float64,
                     multi_relatio::Float64)
    C = Matrix{Float64}(undef,skus,skus) .= 0.0
    D = Matrix{Int64}(undef,skus,skus)   .= 0
    if max_dependence > 0.0
        for i = 1:ceil(skus*ind_chance)
            ind = rand(1:skus)
            D[ind,ind] = 1
        end
        finished = 0
        while finished == 0
            group_size = rand(1:max_group_size)
            test_correlation = transpose(sum(D,dims=1)) .+ sum(D,dims=2)
            free_skus = findall(==(0), test_correlation)
            if group_size <= length(free_skus)
                members = Vector{Int64}(undef,group_size)
                for i = 1:group_size
                    members[i] = rand(1:length(free_skus))
                end
                for i = 1:group_size
                    for j = 1:group_size
                        if rand() < multi_relatio 
                            if i != j
                                strength = rand(Uniform(min_dependence,max_dependence),1)[1]
                                local_i = free_skus[members[i]][1]
                                local_j = free_skus[members[j]][1]
                                D[local_i,local_j] = 1
                                C[local_i,local_j] = strength
                                if rand() > one_direction
                                    D[local_j,local_i] = 1
                                    C[local_j,local_i] = strength
                                end
                            end
                        end
                    end
                end
            else
                finished = 1
            end
        end
        for i in randperm(skus)[1:ceil(Int64,skus*group_link)]
            j = rand(1:skus)
            strength = rand(Uniform(min_dependence,max_dependence),1)[1]
            C[i,j] = strength
            D[i,j] = 1
            if rand() > one_direction
                C[j,i] = strength
                D[j,i] = 1
            end
        end
    end
    C = dropzeros(sparse(C))
    D = dropzeros(sparse(D))
    trans = spzeros(0,skus)
    if rem(orders,1000) != 0 || orders <= 10000
        divide = 1
    else
        divide = orders/10000
    end
    for part = 1:divide
        transactions = spzeros(round(Int64,orders/divide),skus)
        for i = 1:round(Int64,orders/divide)
            already_allocated = 0
            skus_order = 1 + floor(abs(rand(Normal(0,3))))
            while already_allocated < skus_order 
                new_sku = skus + 1
                while new_sku > skus
                    new_sku = 1+ floor(Int64, abs(rand(Normal(0,skus/2.5))))
                end
                already_allocated += 1
                transactions[i,new_sku] = 1
                if rand() > ind_chance
                    filled_D = findnz(D[:,new_sku])[1]
                    if length(filled_D) > 0
                        for j in randperm(length(filled_D))
                        #for j in randperm(skus)
                            #if D[new_sku,j] == 1
                            if transactions[i,filled_D[j]] == 0
                                if already_allocated < skus_order
                                    if rand() < C[filled_D[j],new_sku]
                                        transactions[i,filled_D[j]] = 1
                                        already_allocated += 1
                                    end
                                else
                                    break
                                end
                            end
                            #end
                        end
                    end
                end
            end
        end
        #transactions = sparse(transactions')
        trans = vcat(trans,transactions)
    end
    return trans
end
