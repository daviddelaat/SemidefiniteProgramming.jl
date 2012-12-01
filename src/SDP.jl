module SDP

import Base.print, Base.show

export Rel, SparseSDP, sdpa, sdpa_qd, sdpa_gmp, csdp, print, show

type Rel{T<:Number}
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
    println("Sparse SDP:")
    print(io, sdp)
end

function sdpa(sdp::SparseSDP, app)
    datafname, dataio = mktemp()
    print(dataio, sdp)
    flush(dataio)
    for l in each_line(`$app -ds $datafname -o /dev/null`)
        if begins_with(l, "objValPrimal = ")
            close(dataio)
            return float(split(l, " = ")[2])
        end      
    end
    close(dataio)
    return nothing
end

sdpa(sdp::SparseSDP) = sdpa(sdp, "sdpa")

sdpa_qd(sdp::SparseSDP) = sdpa(sdp, "sdpa_qd")

sdpa_gmp(sdp::SparseSDP) = sdpa(sdp, "sdpa_gmp")

function csdp(sdp::SparseSDP)
    datafname, dataio = mktemp()
    print(dataio, sdp)
    flush(dataio)    
    for l in each_line(`csdp $datafname /dev/null`)
        if begins_with(l, "Primal objective value: ")
            close(dataio)
            return float(split(l, ": ")[2])
        end
    end
    close(dataio)
    return nothing
end

end
