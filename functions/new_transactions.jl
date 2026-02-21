# RANDOMTRANS: function that generates random transactions with dependencies between
#              SKUs based on the values of max_groupsize and max_dependence
function RANDOMTRANS(skus::Int64,
                     orders::Int64,
                     skus_in_order::Float64,
                     sku_frequency::Float64,
                     max_group_size::Int64,
                     min_dependence::Float64,
                     max_dependence::Float64,
                     group_link::Float64,
                     ind_chance::Float64,
                     one_direction::Float64,
                     multi_relatio::Float64)

    C_rows = Int64[]
    C_cols = Int64[]
    C_vals = Float64[]
    taken = zeros(Bool, skus)

    if max_dependence > 0.0
        for i = 1:ceil(Int64, skus * ind_chance)
            ind = rand(1:skus)
            push!(C_rows, ind); push!(C_cols, ind); push!(C_vals, 1.0)
            taken[ind] = true
        end
        finished = 0
        while finished == 0
            group_size = rand(1:max_group_size)
            free_skus = findall(.!taken)
            if group_size <= length(free_skus)
                members = Vector{Int64}(undef, group_size)
                for i = 1:group_size
                    members[i] = rand(1:length(free_skus))
                end
                for i = 1:group_size
                    for j = 1:group_size
                        if rand() < multi_relatio
                            if i != j
                                strength = rand(Uniform(min_dependence, max_dependence), 1)[1]
                                local_i = free_skus[members[i]]
                                local_j = free_skus[members[j]]
                                push!(C_rows, local_i); push!(C_cols, local_j); push!(C_vals, strength)
                                taken[local_i] = true
                                taken[local_j] = true
                                if rand() > one_direction
                                    push!(C_rows, local_j); push!(C_cols, local_i); push!(C_vals, strength)
                                end
                            end
                        end
                    end
                end
            else
                finished = 1
            end
        end
        for i in randperm(skus)[1:ceil(Int64, skus * group_link)]
            j = rand(1:skus)
            strength = rand(Uniform(min_dependence, max_dependence), 1)[1]
            push!(C_rows, i); push!(C_cols, j); push!(C_vals, strength)
            if rand() > one_direction
                push!(C_rows, j); push!(C_cols, i); push!(C_vals, strength)
            end
        end
    end

    C = sparse(C_rows, C_cols, C_vals, skus, skus)
    C_rows = nothing
    C_cols = nothing
    C_vals = nothing
    taken = nothing

    trans = spzeros(Bool, 0, skus)
    if rem(orders, 10000) != 0 || orders <= 10000
        divide = 1
    else
        divide = orders / 10000
    end
    for part = 1:divide
        transactions = spzeros(Bool, round(Int64, orders / divide), skus)
        for i = 1:round(Int64, orders / divide)
            already_allocated = 0
            skus_order = 1 + rand(Geometric(1 / skus_in_order))
            while already_allocated < skus_order
                if sku_frequency == 0
                    new_sku = rand(1:skus)
                else
                    new_sku = skus + 1
                    while new_sku > skus
                        new_sku = 1 + floor(Int64, abs(rand(Normal(0, skus / sku_frequency))))
                    end
                end
                already_allocated += 1
                transactions[i, new_sku] = 1
                if rand() > ind_chance
                    filled_D = findnz(C[:, new_sku])[1]
                    if length(filled_D) > 0
                        for j in randperm(length(filled_D))
                            if transactions[i, filled_D[j]] == 0
                                if already_allocated < skus_order
                                    if rand() < C[filled_D[j], new_sku]
                                        transactions[i, filled_D[j]] = 1
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
        trans = vcat(trans, transactions)
    end
    return trans
end
