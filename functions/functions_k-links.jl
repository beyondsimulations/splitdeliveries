function INL(j::Int64,
            w1::Int64,
            L::Array{Float64,2},
            X::Array{Bool,2})
    out = 0.0
    for q = 1:size(L,1)
        if  X[q,w1] == 1 &&  X[j,w1] == 1 && q != j
            out += L[j,q]
        end
    end
    return out::Float64
end

function OUTL(j::Int64,
              w1::Int64,
              w2::Int64,
              L::Array{Float64,2},
              X::Array{Bool,2})
    out = 0.0
    for q = 1:size(L,1)
        if  X[q,w2] == 1 &&  X[j,w1] == 1 && q != j
            out += L[j,q]
        end
    end
    return out::Float64
end 

function DL(j::Int64,
            w1::Int64,
            w2::Int64,
            L::Array{Float64,2},
            X::Array{Bool,2})
    OUTL(j,w1,w2,L,X) - INL(j,w1,L,X)
end

function PAYOFF(j::Int64,
                q::Int64,
                w1::Int64,
                w2::Int64,
                L::Array{Float64,2},
                X::Array{Bool,2})
    DL(j,w1,w2,L,X) + DL(q,w2,w1,L,X) - 2 * L[j,q]
end

function LINKS(trans::SparseMatrixCSC{Bool, Int64},
               ov::Vector{Float64})
    L = zeros(Float64,size(trans,2),size(trans,2))
    for j = 2:size(trans,2)
        for q = 1:j-1
            L[j,q] = sum(@view(trans[:,j]) .* @view(trans[:,q]) .* @view(ov[:]))
            L[q,j] = L[j,q]
        end
    end
    return L::Matrix{Float64}
end

function LINKADJUST(trans::SparseMatrixCSC{Bool, Int64})
    ov = zeros(Float64,size(trans,1))
    transposed = trans'
    for i = 1:size(trans,1)
        if sum(transposed[:,i]) > 1
            ov[i] = sum(transposed[:,i])-1/(factorial(big(sum(transposed[:,i])))/factorial(2)*factorial(big(sum(transposed[:,i])-2)))
        end
    end
    return ov::Vector{Float64}
end

function LWbest(L::Matrix{Float64})
    LW_best = 0
    for j = 1:size(L,2)
        for q = j+1:size(L,2)
            LW_best += L[j,q]
        end
    end
    return LW_best::Int64
end

function LW(L::Matrix{Float64},
            X::Array{Bool,2})
    out = 0
    for k = 1:(size(X,2))
        for j = 1:size(X,1)
            for q = j+1:size(X,1)
                out += L[j,q] * X[j,k] * X[q,k]
            end
        end
    end
    return out::Float64
end

function STRATEGY1!(X::Matrix{Bool},
                    m::Int64,
                    L::Matrix{Float64},
                    capacity::Vector{Int64},
                    stop::Int64)
    dl_arr = zeros(Float64,size(X,1),size(X,2))
    for i in 1:size(X,1)
        if X[i,m] == 1
            for g = 1:size(X,2)
                if sum(X[:,g]) < capacity[g] && m!=g
                    dl_arr[i,g] = DL(i,m,g,L,X)
                end
            end
        end
    end
    if findmax(dl_arr)[1] > 0
        X[getindex(findmax(dl_arr)[2],1),m] = 0
        X[getindex(findmax(dl_arr)[2],1),getindex(findmax(dl_arr)[2],2)] = 1
        stop = 0
    end
end

function STRATEGY2!(X::Matrix{Bool},
                    m::Int64,
                    L::Matrix{Float64},
                    capacity::Vector{Int64},
                    stop::Int64)
    pay_arr = zeros(Float64,size(X,1),size(X,1),size(capacity,1))
    for i in 1:size(X,1)
        if X[i,m] == 1 
            for j in 1:size(X,1)
                if X[j,m] == 0
                    for g = 1:size(X,2)
                        if X[j,g] == 1 && X[i,g] == 0
                            pay_arr[i,j,g] = PAYOFF(i,j,m,g,L,X)
                        end
                    end
                end
            end
        end
    end
    if findmax(pay_arr)[1] > 0
        X[getindex(findmax(pay_arr)[2],1),m] = 0
        X[getindex(findmax(pay_arr)[2],1),getindex(findmax(pay_arr)[2],3)] = 1
        X[getindex(findmax(pay_arr)[2],2),getindex(findmax(pay_arr)[2],3)] = 0
        X[getindex(findmax(pay_arr)[2],2),m] = 1
        stop = 0
    end
end