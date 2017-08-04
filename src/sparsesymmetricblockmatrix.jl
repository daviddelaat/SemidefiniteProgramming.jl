mutable struct SparseSymmetricBlockMatrix{T<:Number}
    blocks::Dict{Any,SparseSymmetricMatrix{T}}
end

SparseSymmetricBlockMatrix(T::Type) = SparseSymmetricBlockMatrix(Dict{Any,SparseSymmetricMatrix{T}}())

SparseSymmetricBlockMatrix() = SparseSymmetricBlockMatrix(Float64)

blocks(m::SparseSymmetricBlockMatrix) = m.blocks

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, v::T, bi, i, j)
    if haskey(blocks(bm), bi)
        blocks(bm)[bi][i, j] = v
    else
        m = SparseSymmetricMatrix(T)
        m[i, j] = v
        blocks(bm)[bi] = m
    end
end

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, m::SparseSymmetricMatrix{T}, bi)
    blocks(bm)[bi] = m
end

function getindex{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, bi, i, j)
    get(blocks(bm), bi, SparseSymmetricMatrix(T))[i, j]
end

getindex{T}(bm::SparseSymmetricBlockMatrix{T}, bi) = get(blocks(bm), bi, SparseSymmetricMatrix(T))

start(m::SparseSymmetricBlockMatrix) = start(blocks(m))
next(m::SparseSymmetricBlockMatrix, state) = next(blocks(m), state)
done(m::SparseSymmetricBlockMatrix, state) = done(blocks(m), state)

copy(m::SparseSymmetricBlockMatrix) = SparseSymmetricBlockMatrix(copy(blocks(m)))

function *{T}(x::T, m::SparseSymmetricBlockMatrix{T})
    p = SparseSymmetricBlockMatrix(T)
    for k in keys(blocks(m))
        p[k] = x * m[k]
    end
    p
end

function show{T<:Number}(io::IO, m::SparseSymmetricBlockMatrix{T})
    println(io, typeof(m))
    for (blockindex, matrix) in blocks(m)
        for ((i, j), v) in matrix.entries
            println(io, "$blockindex\t$i\t$j\t$v")
        end
    end
end

==(A::SparseSymmetricBlockMatrix, B::SparseSymmetricBlockMatrix) = blocks(A) == blocks(B)

function isapprox(A::SparseSymmetricBlockMatrix, B::SparseSymmetricBlockMatrix)
    for (bi, bm) in A
        if !isapprox(bm, B[bi])
            return false
        end
    end
    for (bi, bm) in B
        if !isapprox(bm, A[bi])
            return false
        end
    end
    true
end
