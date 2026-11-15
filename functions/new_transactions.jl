# RANDOMTRANS: Generates synthetic transactional datasets with controlled
#              dependency structures for benchmarking warehouse allocation heuristics.

# Largest and mean group size for a catalog of skus SKUs: (max_group_size, mean_group_size)
function GROUPSIZES(
    skus::Int64,
    group_size_scaling::Real,
    group_size_cap::Real,
    group_size_min::Real,
    mean_group_divisor::Real,
    mean_group_min::Real,
)
    max_group_size = ceil(
        Int64, max(min(skus / group_size_scaling, group_size_cap), group_size_min)
    )
    mean_group_size = ceil(Int64, max(max_group_size / mean_group_divisor, mean_group_min))
    return max_group_size, mean_group_size
end

function RANDOMTRANS(
    skus::Int64,
    orders::Int64,
    mean_order_size::Float64,
    min_order_size::Int64,
    nbd_dispersion::Float64,
    sku_frequency_mode::Symbol,
    zipf_exponent::Float64,
    max_group_size::Int64,
    mean_group_size::Int64,
    ratio_strong::Float64,
    ratio_medium::Float64,
    dep_strength_strong::Tuple{Float64,Float64},
    dep_strength_medium::Tuple{Float64,Float64},
    group_link::Float64,
    one_direction::Float64,
    multi_relatio::Float64,
    dep_activation_prob::Float64,
)

    # =========================================================================
    # Phase 1: Build dependency structure
    # =========================================================================

    C_rows = Int64[]
    C_cols = Int64[]
    C_vals = Float64[]
    group_sizes = Int64[]

    n_strong = floor(Int64, skus * ratio_strong)
    n_medium = floor(Int64, skus * ratio_medium)

    # Step 1: Partition SKUs into pools
    pool = fill(:independent, skus)
    shuffled_indices = randperm(skus)
    for i in 1:n_strong
        pool[shuffled_indices[i]] = :strong
    end
    for i in (n_strong + 1):(n_strong + n_medium)
        pool[shuffled_indices[i]] = :medium
    end

    # Step 2: Build groups within each dependency pool
    if ratio_strong > 0.0 || ratio_medium > 0.0
        for (tier, dep_strength) in
            [(:strong, dep_strength_strong), (:medium, dep_strength_medium)]
            remaining = findall(pool .== tier)
            if isempty(remaining)
                continue
            end

            while length(remaining) >= 2
                # Group size from geometric distribution (many small, few large)
                gsize = min(
                    2 + rand(Geometric(1.0 / (mean_group_size - 1))),
                    max_group_size,
                    length(remaining),
                )

                # Sample group members without replacement (using randperm)
                perm = randperm(length(remaining))
                member_indices = sort(perm[1:gsize])
                members = remaining[member_indices]

                # Remove members from remaining pool (reverse order to preserve indices)
                deleteat!(remaining, member_indices)

                push!(group_sizes, gsize)

                # Designate first member as anchor/hub
                anchor = members[1]
                peripherals = members[2:end]

                # Anchor connects to ALL peripherals
                for p in peripherals
                    strength = rand(Uniform(dep_strength[1], dep_strength[2]))
                    push!(C_rows, p);
                    push!(C_cols, anchor);
                    push!(C_vals, strength)
                    if rand() > one_direction
                        push!(C_rows, anchor);
                        push!(C_cols, p);
                        push!(C_vals, strength)
                    end
                end

                # Peripheral-to-peripheral connections (sparser)
                for i in 1:length(peripherals)
                    for j in (i + 1):length(peripherals)
                        if rand() < multi_relatio
                            strength = rand(Uniform(dep_strength[1], dep_strength[2]))
                            push!(C_rows, peripherals[j]);
                            push!(C_cols, peripherals[i]);
                            push!(C_vals, strength)
                            if rand() > one_direction
                                push!(C_rows, peripherals[i]);
                                push!(C_cols, peripherals[j]);
                                push!(C_vals, strength)
                            end
                        end
                    end
                end
            end
        end

        # Step 3: Cross-group bridge links
        n_bridges = ceil(Int64, skus * group_link)
        candidates = findall(pool .!= :independent)
        if length(candidates) >= 2
            min_strength = min(dep_strength_strong[1], dep_strength_medium[1])
            max_strength = max(dep_strength_strong[2], dep_strength_medium[2])
            for _ in 1:n_bridges
                i = rand(candidates)
                j = rand(candidates)
                if i != j
                    strength = rand(Uniform(min_strength, max_strength))
                    push!(C_rows, i);
                    push!(C_cols, j);
                    push!(C_vals, strength)
                    if rand() > one_direction
                        push!(C_rows, j);
                        push!(C_cols, i);
                        push!(C_vals, strength)
                    end
                end
            end
        end
    end

    C = sparse(C_rows, C_cols, C_vals, skus, skus)
    # Clamp values to [0,1] — sparse() sums duplicates from overlapping bridge/group edges
    C.nzval .= min.(C.nzval, 1.0)
    C_rows = nothing
    C_cols = nothing
    C_vals = nothing

    # =========================================================================
    # Phase 2: Generate transactions
    # =========================================================================

    # Step 1: Pre-compute SKU selection distribution
    if sku_frequency_mode == :zipf
        rank_to_sku = randperm(skus)
        raw_weights = [1.0 / k^zipf_exponent for k in 1:skus]
        sku_weights = zeros(Float64, skus)
        total_weight = sum(raw_weights)
        for k in 1:skus
            sku_weights[rank_to_sku[k]] = raw_weights[k] / total_weight
        end
        # Cumulative selection probabilities in SKU order, summed one by one as
        # rand(Categorical(sku_weights)) does: the seed SKU for a uniform draw u is the
        # first SKU whose cumulative probability exceeds u, or the last SKU
        sku_cumulative = cumsum(sku_weights)
    else
        sku_cumulative = nothing
    end

    # Step 2: Pre-compute order size distribution (Negative Binomial)
    effective_mean = mean_order_size - min_order_size
    if effective_mean <= 0.0
        effective_mean = 0.01
    end
    p_nbd = nbd_dispersion / (nbd_dispersion + effective_mean)
    order_size_dist = NegativeBinomial(nbd_dispersion, p_nbd)

    # Step 3: Generate orders; the SKUs of order i are
    # order_skus[order_start[i]:(order_start[i + 1] - 1)] in ascending order
    order_skus = Int64[]
    order_start = Vector{Int64}(undef, orders + 1)
    order_start[1] = 1

    for i in 1:orders
        # Draw order size from NBD, capped at total number of SKUs
        order_size = min(min_order_size + rand(order_size_dist), skus)
        items_in_order = Set{Int64}()

        while length(items_in_order) < order_size
            # Draw seed SKU from frequency distribution
            if sku_frequency_mode == :zipf
                seed = min(searchsortedlast(sku_cumulative, rand()) + 1, skus)
            else
                seed = rand(1:skus)
            end

            push!(items_in_order, seed)

            # Multi-hop dependency propagation via BFS queue
            # Each newly added SKU can trigger its own dependents
            queue = [seed]
            while !isempty(queue) && length(items_in_order) < order_size
                current = popfirst!(queue)
                dep_range = nzrange(C, current)
                dep_skus = view(rowvals(C), dep_range)
                dep_strengths = view(nonzeros(C), dep_range)

                if length(dep_skus) > 0 && rand() < dep_activation_prob
                    for j in randperm(length(dep_skus))
                        if length(items_in_order) >= order_size
                            break
                        end
                        if dep_skus[j] ∉ items_in_order && rand() < dep_strengths[j]
                            push!(items_in_order, dep_skus[j])
                            push!(queue, dep_skus[j])
                        end
                    end
                end
            end
        end

        # Write the SKUs of the order in ascending order
        append!(order_skus, items_in_order)
        sort!(view(order_skus, order_start[i]:length(order_skus)))
        order_start[i + 1] = length(order_skus) + 1
    end

    # One column per order, turned into one row per order
    trans = permutedims(
        SparseMatrixCSC(skus, orders, order_start, order_skus, fill(true, length(order_skus)))
    )

    return trans, C, group_sizes
end
