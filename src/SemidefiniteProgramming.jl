module SemidefiniteProgramming

import 
    Base.print, 
    Base.show,
    Base.getindex,
    Base.setindex!
    
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
    parsesdpasparse,
    SDPSolver,
    SDPA,
    SDPAQD,
    SDPAGMP,
    CSDP,
    solve

   
    
    
type SparseSymmetricMatrix{T<:Number}
    m::Dict{(Int, Int), T}
end

SparseSymmetricMatrix(T) = SparseSymmetricMatrix{T}(Dict{(Int,Int), T)
SparseSymmetricMatrix() = SparseSymmetricMatrix(Float64)

function setindex!{T<:Number}(m::SparseSymmetricMatrix{T}, v::T, i::Int, j::Int)
    if i <= j
        m[(i, j)] = v
    else
        m[(j, i)] = v
    end
end

function getindex{T<:Number}(m::SparseSymmetricMatrix{T}, i::Int, j::Int)
    if i <= j
        get(m, (i, j), zero(T))
    else
        get(m, (j, i), zero(T))
    end
end




type SparseSymmetricBlockMatrix{T<:Number}
    blocks::Dict{Int, SparseSymmetricMatrix{T}}
end

SparseSymmetricBlockMatrix(T::Number) = SparseSymmetricBlockMatrix(Dict{Int, SparseSymmetricMatrix{T}}())
SparseSymmetricBlockMatrix() = SparseSymmetricBlockMatrix(Float64)

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, v::T, bi::Int, i::Int, j::Int)
    get(bm.blocks, bi, SparseSymmetricMatrix(T))[i, j] = v
end

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, m::SparseSymmetricMatrix{T}, bi::Int)
    bm.blocks[bi] = m
end

function getindex{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, bi::Int, i::Int, j::Int)
    get(bm.blocks, bi, SparseSymmetricMatrix(T))[i, j]
end




type SparseSDP{T<:Number}
    c::SparseSymmetricBlockMatrix{T}
    as::Dict{Int, SparseSymmetricBlockMatrix{T}}
    bs::Dict{Int, T}
end

SparseSDP(T::Number) = SparseSDP(SparseSymmetricBlockMatrix(T), 
                                 Dict{Int, SparseSymmetricBlockMatrix{T}}(),
                                 Dict{Int, T}())

SparseSDP() = SparseSDP(Float64)

setc!{T<:Number}(sdp::SparseSDP{T}, c::SparseBlockMatrix{T}) = sdp.c = c

function setc!{T<:Number}(sdp::SparseSDP{T}, bi::Int, m::SparseSymmetricMatrix{T})
    sdp.c[bi] = m
end

function setc!{T<:Number}(sdp::SparseSDP{T}, bi::Int, i::Int, j::Int, v::T)
    sdp.c[bi, i, j] = v
end

getc(sdp::SparseSDP) = sdp.c

function seta!{T<:Number}(sdp::SparseSDP{T}, ri::Int, bi::Int, i::Int, j::Int, v::T)
    if has(sdp.as, ri)
        sdp.as[ri][bi, i, j] = v
    else 
    sdp.tcs[ri] = tc
end

geta(sdp::SparseSDP, ri::Int) = sdp.tcs[ri]

geta(sdp::SparseSDP, ri::Int, bi::Int) = sdp. 




function blocksizes(sdp::SparseSDP)
    bsizes = Dict{Int, Int}()
    mbi = 0
    for (ri, bmat) in sdp.as
        for 
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

function m(sdp::SparseSDP)
    mm = 0
    for ri in keys(sdp.as)
        if ri > mm
            mm = ri
        end
    end
    for ri in keys(sdp.bs)
        if ri > mm
            mm = ri
        end
    end
    mm
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

function parsesdpasparse(io::IO)
    sdp = SparseSDP()
    
    s = 1
    for l in readlines(io)
        println(l)
        if !beginswith(l, '#')
            if s == 4
                t = split(l) 
                for i = 1:length(t)
                     sdp[i] = parsefloat(t[i])
                end
            elseif s > 4
                t = split(l)
                sdp[parseint(t[1]), parseint(t[2]), parseint(t[3]), parseint(t[4])] = parsefloat(t[5])
            end 
            s += 1
        end
    end
    
    sdp
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
