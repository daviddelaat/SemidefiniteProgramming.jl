mutable struct SparseSDPSolution{T<:Number}
    obj::T
    dualobj::T
    primalmatrix::SparseSymmetricBlockMatrix{T}
    dualvector::Dict{Any,T}
    dualmatrix::SparseSymmetricBlockMatrix{T}
end

SparseSDPSolution(T) = SparseSDPSolution(zero(T), zero(T), Dict(), Dict(), Dict())

SparseSDPSolution() = SparseSDPSolution(Float64)

SparseSDPSolution{T<:Number}(obj::T, dualobj::T) = SparseSDPSolution{T}(obj, 
                                                                        dualobj, 
                                                                        SparseSymmetricBlockMatrix(T), 
                                                                        Dict{Any,T}(), 
                                                                        SparseSymmetricBlockMatrix(T))

setobj!{T<:Number}(sol::SparseSDPSolution{T}, v::T) = sol.obj = v

obj(sol::SparseSDPSolution) = sol.obj

setdualobj!{T<:Number}(sol::SparseSDPSolution{T}, v::T) = sol.dualobj = v

dualobj(sol::SparseSDPSolution) = sol.dualobj

primalmatrix(sol::SparseSDPSolution) = sol.primalmatrix

dualvector(sol::SparseSDPSolution) = sol.dualvector

dualmatrix(sol::SparseSDPSolution) = sol.dualmatrix

function invdict(d::Dict)
    r = Dict()
    for (k, v) in d
        r[v] = k
    end
    r
end

function readcsdpoutput!(io::IO, sol::SparseSDPSolution, cm::Dict, bm::Dict, ems::Dict)
    icm = invdict(cm)
    ibm = invdict(bm)
    iems = Dict()
    for (k, v) in ems
        iems[ibm[k]] = invdict(v)
    end
    
    l = readline(io)
    ls = split(l)
    for i = 1:length(ls)
        dualvector(sol)[icm[i]] = float(ls[i])
    end
    
    for l in readlines(io)
        ls = split(l)
        onetwo = parse(Int, ls[1])
        blockindex = ibm[parse(Int, ls[2])]
        rowindex = iems[blockindex][parse(Int, ls[3])]
        columnindex = iems[blockindex][parse(Int, ls[4])]
        value = parse(Float64, ls[5])
        if onetwo == 2
            block = get(primalmatrix(sol).blocks, blockindex, SparseSymmetricMatrix())
            block[rowindex, columnindex] = value
            primalmatrix(sol).blocks[blockindex] = block
        else
            block = get(dualmatrix(sol).blocks, blockindex, SparseSymmetricMatrix())
            block[rowindex, columnindex] = value
            dualmatrix(sol).blocks[blockindex] = block
        end
    end
    
    sol
end
