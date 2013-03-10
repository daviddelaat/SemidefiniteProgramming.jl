module SemidefiniteProgramming

import 
    Base.print, 
    Base.show

export 
    Rel, 
    SparseSDP, 
    SDPSolver,
    SDPA,
    SDPAQD,
    SDPAGMP,
    CSDP,
    solve


    
    
    
# Sparse SDPA format
    
immutable Rel{T<:Number}
    nr::Int
    blocknr::Int
    i::Int
    j::Int
    val::T
end

type SparseSDP{T<:Number}
    blocksizes::Vector{Int}
    rhs::Vector{T}
    rels::Vector{Rel{T}}
end

function SparseSDP{T<:Number}(rhs::Vector{T}, rels::Vector{Rel{T}})
    m = 0
    for rel in rels
        if rel.blocknr > m
            m = rel.blocknr
        end
    end
    blocksizes = zeros(Int, m)
    for rel in rels
        m = max(rel.i, rel.j)
        if m > blocksizes[rel.blocknr]
            blocksizes[rel.blocknr] = m
        end
    end
    SparseSDP(blocksizes, rhs, rels)
end

function print(io::IO, sdp::SparseSDP)
    println(io, length(sdp.rhs))
    println(io, length(sdp.blocksizes))
    
    for b in sdp.blocksizes[1:end-1]
        print(io, b)
        print(io, " ")
    end
    println(io, sdp.blocksizes[end])
    
    for r in sdp.rhs[1:end-1]
        print(io, r)
        print(io, " ")
    end
    println(io, sdp.rhs[end])
    
    for rel in sdp.rels
        print(io, rel.nr)
        print(io, " ")
        print(io, rel.blocknr)
        print(io, " ")
        print(io, rel.i)
        print(io, " ")
        print(io, rel.j)
        print(io, " ")
        print(io, rel.val)
        println(io)  
    end
end

function show(io::IO, sdp::SparseSDP)
    println("$(length(sdp.blocksizes))-block sparse semidefinite program:")
    print(io, sdp)
end






# Solvers

abstract SDPSolver

abstract SDPAGEN <: SDPSolver

type SDPA <: SDPAGEN
    speed::Int
    executable::ASCIIString
end
SDPA{T<:Integer}(speed::T) = SDPA(speed, "sdpa")
SDPA() = SDPA(0)

type SDPAQD <: SDPAGEN
    speed::Int
    executable::ASCIIString
end
SDPAQD{T<:Integer}(speed::T) = SDPAQD(speed, "sdpa_qd")
SDPAQD() = SDPAQD(0)

type SDPAGMP <: SDPAGEN
    speed::Int
    executable::ASCIIString
end
SDPAGMP(speed::Int) = SDPAGMP(speed, "sdpa_gmp")
SDPAGMP() = SDPAGMP(0)

type CSDP <: SDPSolver
    executable::ASCIIString
end
CSDP() = CSDP("csdp")

function solve{T<:SDPAGEN}(sdp::SparseSDP, solver::T)
    datafname, dataio = mktemp()
    print(dataio, sdp)
    flush(dataio)
    for l in each_line(`$(solver.executable) -ds $datafname -o /dev/null`)
        if begins_with(l, "objValPrimal = ")
            close(dataio)
            return float(split(l, " = ")[2])
        end      
    end
    close(dataio)
    return nothing
end

function solve(sdp::SparseSDP, solver::CSDP)
    datafname, dataio = mktemp()
    print(dataio, sdp)
    flush(dataio)    
    for l in each_line(`$(solver.executable) $datafname /dev/null`)
        if begins_with(l, "Primal objective value: ")
            close(dataio)
            return float(split(l, ": ")[2])
        end
    end
    close(dataio)
    return nothing
end

end
