## import packages
include("load_packages.jl")
using SparseArrays

T = RANDOMTRANS(1000,100000,100,100,0.0,0.0)

function COAPPEARENCE(trans::Matrix{Int64})
    Q = Octavian.matmul(transpose(trans),trans)
    for i = 1:size(Q,1)
        Q[i,i] = 0
    end
    return Q::Array{Int64,2}
end

#function COAPPEARENCE(trans::SparseMatrixCSC{Int64,Int64})#
#    Q = zeros(Int64,size(trans,2),size(trans,2))
#    for i = 2:size(Q,1)
#        for j =1:i
#            for k = 1:size(trans,2)
#                if trans[k,i] & trans[k,j] == 1
#                    Q[i,j] += 1
#                end
#            end
#        end
#    end
#    return Q::Array{Int64,2}
#end

function COAPPEARENCE(trans::SparseMatrixCSC{Int64,Int64})
    t_sparse = dropzeros(sparse(trans))
    Q = t_sparse'*t_sparse
    for i = 1:size(Q,1)
        Q[i,i] = 0
    end
    Q = Matrix(Q)
    return Q
end

@time q = COAPPEARENCE(T)
@time q = COAPPEARENCE(T_sparse)