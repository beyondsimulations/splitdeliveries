function INL(j::Int64, w1::Int64, L::Array{Float64,2}, X::Array{Bool,2})
    out = 0.0
    @inbounds for q in axes(L, 1)
        if X[q, w1] & X[j, w1] & (q != j)
            out += L[j, q]
        end
    end
    return out::Float64
end

function OUTL(j::Int64, w1::Int64, w2::Int64, L::Array{Float64,2}, X::Array{Bool,2})
    out = 0.0
    @inbounds for q in axes(L, 1)
        if X[q, w2] & X[j, w1] & (q != j)
            out += L[j, q]
        end
    end
    return out::Float64
end

function DL(j::Int64, w1::Int64, w2::Int64, L::Array{Float64,2}, X::Array{Bool,2})
    OUTL(j, w1, w2, L, X) - INL(j, w1, L, X)
end

function PAYOFF(
    j::Int64, q::Int64, w1::Int64, w2::Int64, L::Array{Float64,2}, X::Array{Bool,2}
)
    DL(j, w1, w2, L, X) + DL(q, w2, w1, L, X) - 2 * L[j, q]
end

function SPARSELINKS(trans::SparseMatrixCSC{Bool,Int64}, ov::Vector{Float64})
    # L[j, q] = sum over orders m of ov[m] * trans[m, j] * trans[m, q], zero diagonal
    L = trans' * (Diagonal(ov) * trans)
    CLEANPRINCIPLE!(L)
    return L
end

function LINKS(trans::SparseMatrixCSC{Bool,Int64}, ov::Vector{Float64})
    return Matrix(SPARSELINKS(trans, ov))
end

function LINKADJUST(trans::SparseMatrixCSC{Bool,Int64})
    ov = zeros(Float64, size(trans, 1))
    order_size = vec(sum(trans; dims = 2))
    for i in axes(trans, 1)
        n = order_size[i]
        if n > 1
            # contribution effect coefficient of Zhu et al. (2021): an order with n
            # SKUs can cause at most n-1 splits, spread over binomial(n,2) pairs
            ov[i] = (n - 1) / binomial(n, 2)
        end
    end
    return ov::Vector{Float64}
end

function LWbest(L::Matrix{Float64})
    LW_best = 0
    for j in axes(L, 2)
        for q in (j + 1):size(L, 2)
            LW_best += L[j, q]
        end
    end
    return LW_best::Int64
end

function LW(L::Matrix{Float64}, X::Array{Bool,2})
    out = 0
    for k in 1:(size(X, 2))
        for j in axes(X, 1)
            for q in (j + 1):size(X, 1)
                out += L[j, q] * X[j, k] * X[q, k]
            end
        end
    end
    return out::Float64
end

# Relocate SKU i to warehouse g and update S = L * X, wh and load
function RELOCATE!(
    X::Matrix{Bool},
    S::Matrix{Float64},
    wh::Vector{Int64},
    load::Vector{Int64},
    L::Matrix{Float64},
    i::Int64,
    g::Int64,
)
    S[:, wh[i]] .-= view(L, :, i)
    S[:, g] .+= view(L, :, i)
    X[i, wh[i]] = false
    X[i, g] = true
    load[wh[i]] -= 1
    load[g] += 1
    wh[i] = g
    return nothing
end
