type SparseSymmetricMatrix{T<:Number,I}
    entries::Dict{(I, I), T}
    indices::Set{I}
end

SparseSymmetricMatrix(T::Type, I::Type) = SparseSymmetricMatrix{T,I}(Dict{(I,I),T}(), Set{I}())
SparseSymmetricMatrix(T::Type) = SparseSymmetricMatrix(T, Any)
SparseSymmetricMatrix() = SparseSymmetricMatrix(Float64)

entries(m::SparseSymmetricMatrix) = m.entries
indices(m::SparseSymmetricMatrix) = m.indices

function setindex!{T<:Number,I}(m::SparseSymmetricMatrix{T,I}, v::T, i, j)
    if method_exists(isless, (I, I))
        if i <= j
            entries(m)[(i, j)] = v
        else
            entries(m)[(j, i)] = v
        end
    elseif method_exists(hash, (I,)) && method_exists(hash, (I,)) 
        if hash(i) <= hash(j)
            entries(m)[(i, j)] = v
        else
            entries(m)[(j, i)] = v
        end
    else
        entries(m)[(i, j)] = v
        entries(m)[(j, i)] = v
    end
    add!(indices(m), i)
    add!(indices(m), j)
end

function getindex{T<:Number,I}(m::SparseSymmetricMatrix{T,I}, i, j)
    if method_exists(isless, (I, I))
        if i <= j
            get(m.entries, (i, j), zero(T))
        else
            get(m.entries, (j, i), zero(T))
        end
    elseif method_exists(hash, (I,)) && method_exists(hash, (I,)) 
        if hash(i) <= hash(j)
            get(m.entries, (i, j), zero(T))
        else
            get(m.entries, (j, i), zero(T))
        end
    else
        get(m.entries, (i, j), zero(T))
    end
end

start(m::SparseSymmetricMatrix) = start(m.entries)
done(m::SparseSymmetricMatrix, state) = done(m.entries, state)
next(m::SparseSymmetricMatrix, state) = next(m.entries, state)

function size(m::SparseSymmetricMatrix)
    is = Set()
    js = Set()
    for (i, j) in keys(m.entries)
        add!(is, i)
        add!(js, j)
    end
    return max(length(is), length(js))
end
