type SparseSymmetricBlockMatrix{T<:Number,I}
    blocks::Dict{I,SparseSymmetricMatrix{T}}
end

SparseSymmetricBlockMatrix(T::Type, I) = SparseSymmetricBlockMatrix(Dict{I,SparseSymmetricMatrix{T}}())

SparseSymmetricBlockMatrix(T::Type) = SparseSymmetricBlockMatrix(T, Any)

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

function *{T,I}(x::T, m::SparseSymmetricBlockMatrix{T,I})
    p = SparseSymmetricBlockMatrix(T, I)
    for k in keys(blocks(m))
        p[k] = x * m[k]
    end
    p
end
