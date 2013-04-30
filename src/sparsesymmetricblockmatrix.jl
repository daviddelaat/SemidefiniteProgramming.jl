type SparseSymmetricBlockMatrix{T<:Number,I}
    blocks::Dict{I,SparseSymmetricMatrix{T}}
end

SparseSymmetricBlockMatrix(T::Type, I) = SparseSymmetricBlockMatrix(Dict{I,SparseSymmetricMatrix{T}}())
SparseSymmetricBlockMatrix(T::Type) = SparseSymmetricBlockMatrix(T, Any)
SparseSymmetricBlockMatrix() = SparseSymmetricBlockMatrix(Float64)

blocks(m::SparseSymmetricBlockMatrix) = m.blocks

function setindex!{T<:Number,I}(bm::SparseSymmetricBlockMatrix{T,I}, v::T, bi, i, j)
    if haskey(blocks(bm), bi)
        blocks(bm)[bi][i, j] = v
    else
        m = SparseSymmetricMatrix(T)
        m[i, j] = v
        blocks(bm)[bi] = m
    end
end

function setindex!{T<:Number,I}(bm::SparseSymmetricBlockMatrix{T,I}, m::SparseSymmetricMatrix{T}, bi)
    blocks(bm)[bi] = m
end

function getindex(bm::SparseSymmetricBlockMatrix, bi, i, j)
    get(blocks(bm), bi, SparseSymmetricMatrix(T))[i, j]
end

getindex{T}(bm::SparseSymmetricBlockMatrix{T}, bi) = get(blocks(bm), bi, SparseSymmetricMatrix(T))

start(m::SparseSymmetricBlockMatrix) = start(blocks(m))
next(m::SparseSymmetricBlockMatrix, state) = next(blocks(m), state)
done(m::SparseSymmetricBlockMatrix, state) = done(blocks(m), state)
