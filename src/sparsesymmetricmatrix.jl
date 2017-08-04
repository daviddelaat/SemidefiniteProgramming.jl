mutable struct SparseSymmetricMatrix{T<:Number}
    entries::Dict{Tuple{Any, Any}, T}
    indices::Set{Any}
end

SparseSymmetricMatrix(T::Type) = SparseSymmetricMatrix{T}(Dict{Tuple{Any,Any},T}(), Set())

SparseSymmetricMatrix() = SparseSymmetricMatrix(Float64)

entries(m::SparseSymmetricMatrix) = m.entries

indices(m::SparseSymmetricMatrix) = m.indices

function setindex!{T<:Number}(m::SparseSymmetricMatrix{T}, v::T, i, j)
    entries(m)[(i, j)] = v
    push!(indices(m), i)
    push!(indices(m), j)
end

function getindex{T<:Number}(m::SparseSymmetricMatrix{T}, i, j)
    get(m.entries, (i, j), zero(T))
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

copy(m::SparseSymmetricMatrix) = SparseSymmetricMatrix(copy(entries(m)), copy(indices(m)))

function *{T<:Number}(x::T, m::SparseSymmetricMatrix{T})
    p = copy(m)
    for k in keys(entries(m))
        entries(p)[k] *= x
    end
    p
end

function delete!{T<:Number}(m::SparseSymmetricMatrix{T}, i, j)
    delete!(entries(m), (i, j))
    igone = true
    jgone = true
    for (ti, tj) in keys(entries(m))
        if i == ti
            igone = false
        end
        if j == tj
            jgone = false
        end
    end
    if igone
        delete!(indices(m), i)
    end
    if jgone && contains(indices(m), j)
        delete!(indices(m), j)
    end
end

function ==(A::SparseSymmetricMatrix, B::SparseSymmetricMatrix)
    A.entries == B.entries && A.indices == B.indices
end

function isapprox(A::SparseSymmetricMatrix, B::SparseSymmetricMatrix)
    if indices(A) != indices(B)  
        return false
    end
    for ((ri, ci), v) in A
        if !isapprox(v, B[ri, ci]) 
            return false
        end
    end
    for ((ri, ci), v) in B
        if !isapprox(v, A[ri, ci]) 
            return false
        end
    end
    true
end
