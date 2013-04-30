type SparseSDP{T<:Number,I}
    obj::SparseSymmetricBlockMatrix{T}
    cons::Dict{I,SparseSymmetricBlockMatrix{T}}
    rhs::Dict{I,T}
end

SparseSDP(T::Type, I) = SparseSDP(SparseSymmetricBlockMatrix(T), 
                                  Dict{I,SparseSymmetricBlockMatrix{T}}(),
                                  Dict{I,T}())
SparseSDP(T::Type) = SparseSDP(T, Any)
SparseSDP() = SparseSDP(Float64)

setobj!{T<:Number}(sdp::SparseSDP{T}, c::SparseSymmetricBlockMatrix{T}) = obj(sdp) = obj

function setobj!{T<:Number}(sdp::SparseSDP{T}, bi, m::SparseSymmetricMatrix{T})
    obj(sdp)[bi] = m
end

function setobj!{T<:Number}(sdp::SparseSDP{T}, bi, i, j, v::T)
    obj(sdp)[bi, i, j] = v
end

function setobj!{T<:Number}(sdp::SparseSDP{T}, bi, m::Matrix{T})
    n = size(m)[1]
    for i = 1:n
        for j = 1:n
            if m[i, j] != zero(T)
                setobj!(sdp, bi, i, j, m[i, j])
            end
        end
    end
end             

obj(sdp::SparseSDP) = sdp.obj

obj(sdp::SparseSDP, bi) = sdp.obj[bi]

function setcon!{T<:Number}(sdp::SparseSDP{T}, ri, bi, i, j, v::T)
    if haskey(cons(sdp), ri)
        cons(sdp)[ri][bi, i, j] = v
    else
        a = SparseSymmetricBlockMatrix(T)
        a[bi, i, j] = v
        cons(sdp)[ri] = a
    end
end

cons(sdp::SparseSDP) = sdp.cons

getcon(sdp::SparseSDP, ri) = sdp.cons[ri]

setrhs!{T<:Number}(sdp::SparseSDP{T}, ri, v::T) = sdp.rhs[ri] = v

rhs(sdp::SparseSDP) = sdp.rhs

rhs(sdp::SparseSDP, ri) = sdp.rhs[ri]

function rhsdense{T<:Number}(sdp::SparseSDP{T})
    l = zeros(T, ncons(sdp))
    for (ri, v) in rhs(sdp)
        l[ri] = v
    end
    l
end

function ncons(sdp::SparseSDP)
    v = 0
    for ri in keys(cons(sdp))
        if ri > v
            v = ri
        end
    end
    for ri in keys(rhs(sdp))
        if ri > v
            v = ri
        end
    end
    v
end

function blocksizes(sdp::SparseSDP)
    blockdict = Dict{Any,Set}()
    for (bi, block) in obj(sdp)
        indices = get(blockdict, bi, Set())
        union!(indices, block.indices)
        blockdict[bi] = indices
    end
    for a in values(cons(sdp))
        for (bi, block) in a
            indices = get(blockdict, bi, Set())
            union!(indices, block.indices)
            blockdict[bi] = indices
        end
    end
    blocksizedict = Dict{Any,Int}()
    for (bi, entries) in blockdict
        blocksizedict[bi] = length(entries)
    end
    blocksizedict
end
