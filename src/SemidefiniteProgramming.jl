module SemidefiniteProgramming

import 
    Base.print, 
    Base.show
    
export 
    SparseSDP,
    blocksizes,
    setrel!,
    getrel,
    setrhs!,
    getrhs,
    nrels,
    rhs,
    writesdpasparse,
    SDPSolver,
    SDPA,
    SDPAQD,
    SDPAGMP,
    CSDP,
    solve
   
    
    
    

    
    
    
type SparseSDP{T<:Number}
    rhs::Dict{Int, T}
    rels::Dict{(Int, Int, Int, Int), T}
end

SparseSDP{T<:Number}() = SparseSDP(Dict{Int, T}(), Dict{(Int, Int, Int, Int), T}())
SparseSDP() = SparseSDP(Dict{Int, Float64}(), Dict{(Int, Int, Int, Int), Float64}())

setrhs!{T<:Number}(sdp::SparseSDP{T}, ri::Int, v::T) = sdp.rhs[ri] = v       
getrhs(sdp::SparseSDP, ri::Int) = sdp.rhs[ri]

function setrel!{T<:Number}(sdp::SparseSDP{T}, ri::Int, bi::Int, i::Int, j::Int, v::T) 
    if i <= j
        sdp.rels[(ri, bi, i, j)] = v
    else
        sdp.rels[(ri, bi, j, i)] = v
    end
end

setrel!{T<:Number}(sdp::SparseSDP{T}, ri::Int, i::Int, j::Int, v::T) = setrel!(sdp, ri, 1, i, j, v)

function getrel(sdp::SparseSDP, ri::Int, bi::Int, i::Int, j::Int)
    if i <= j
        sdp.rels[(ri, bi, i, j)]
    else
        sdp.rels[(ri, bi, j, i)]
    end
end

getrel(sdp::SparseSDP, ri::Int, i::Int, j::Int) = getrel(sdp, ri, 1, i, j)

function setrel!{T<:Number}(sdp::SparseSDP{T}, ri::Int, bi::Int, m::Matrix{T}) 
    for i = 1:dims(m)[1]
        for j = 1:i
            if m[i, j] != zero(T)
                setrel!(sdp, ri, bi, i, j) = m[i, j]
            end
        end
    end 
end

setrel!{T<:Number}(sdp::SparseSDP{T}, ri::Int, m::Matrix{T}) = setrel!(sdp, ri, 1, m)

function blocksizes(sdp::SparseSDP)
    bsizes = Dict{Int, Int}()
    mbi = 0
    for (ri, bi, i, j) in keys(sdp.rels)
        if bi > mbi
            mbi = bi
        end
        m = max(i, j)
        if m > get(bsizes, bi, 0)
            bsizes[bi] = m
        end
    end
    bsizeslist = zeros(Int, mbi)
    for (bi, bsize) in bsizes
        bsizeslist[bi] = bsize
    end
    bsizeslist
end

function nrels(sdp::SparseSDP)
    c = -1
    for (ri, bi, i, j) in keys(sdp.rels)
        if ri > c
            c = ri
        end
    end
    for ri in keys(sdp.rhs)
        if ri > c
            c = ri
        end
    end
    c    
end

function rhs(sdp::SparseSDP)
    l = zeros(nrels(sdp))
    for (ri, val) in sdp.rhs
        l[ri] = val 
    end
    l
end

function writesdpasparse(io::IO, sdp::SparseSDP)
    bsizes = blocksizes(sdp)    

    println(io, nrels(sdp))
    println(io, length(bsizes))
    
    for b in bsizes[1:end-1]
        print(io, b)
        print(io, " ")
    end
    if length(bsizes) > 0
        println(io, bsizes[end])
    end
    
    rhslist = rhs(sdp)
    for r in rhslist[1:end-1]
        print(io, r)
        print(io, " ")
    end
    if length(rhslist) > 0
        println(io, rhslist[end])
    end
    
    for ((ri, bi, i, j), val) in sdp.rels
        print(io, ri)
        print(io, " ")
        print(io, bi)
        print(io, " ")
        print(io, i)
        print(io, " ")
        print(io, j)
        print(io, " ")
        print(io, val)
        println(io)  
    end
end

function show{T<:Number}(io::IO, sdp::SparseSDP{T})
    bsizes = blocksizes(sdp)
    nr = nrels(sdp)
    
    print(io, "SparseSDP with blocksizes ")
    println(io, bsizes)
    
    for i in 0:nr
        println(io, "Relation $i with right hand side $(get(sdp.rhs, i, zero(T))): ")
        
    end
end







abstract SDPSolver

abstract SDPAGEN <: SDPSolver

immutable SDPA <: SDPAGEN
    speed::Int
    executable::ASCIIString
end
SDPA{T<:Integer}(speed::T) = SDPA(speed, "sdpa")
SDPA() = SDPA(0)

immutable SDPAQD <: SDPAGEN
    speed::Int
    executable::ASCIIString
end
SDPAQD{T<:Integer}(speed::T) = SDPAQD(speed, "sdpa_qd")
SDPAQD() = SDPAQD(0)

immutable SDPAGMP <: SDPAGEN
    speed::Int
    executable::ASCIIString
end
SDPAGMP(speed::Int) = SDPAGMP(speed, "sdpa_gmp")
SDPAGMP() = SDPAGMP(0)

immutable CSDP <: SDPSolver
    executable::ASCIIString
end
CSDP() = CSDP("csdp")

function solve(sdp::SparseSDP, solver::SDPAGEN)
    datafname, dataio = mktemp()
    writesdpasparse(dataio, sdp)
    flush(dataio)
    for l in eachline(`$(solver.executable) -ds $datafname -o /dev/null`)
        if beginswith(l, "objValPrimal = ")
            close(dataio)
            return float(split(l, " = ")[2])
        end      
    end
    close(dataio)
    nothing
end

function solve(sdp::SparseSDP, solver::CSDP)
    datafname, dataio = mktemp()
    writesdpasparse(dataio, sdp)
    flush(dataio)    
    for l in eachline(`$(solver.executable) $datafname /dev/null`)
        if beginswith(l, "Primal objective value: ")
            close(dataio)
            return float(split(l, ": ")[2])
        end
    end
    close(dataio)
    nothing
end

end