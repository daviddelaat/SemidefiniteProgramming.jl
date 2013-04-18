module SemidefiniteProgramming

import 
    Base.print, 
    Base.show,
    Base.getindex,
    Base.setindex!,
    Base.size,
    Base.start,
    Base.next,
    Base.done
    
export 
    SparseSDP,
    blocksizes,
    seta!,
    geta,
    setb!,
    getb,
    setc!,
    getc,
    m,
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
    entries::Dict{(Int, Int), T}
    size::Int
end

SparseSymmetricMatrix(T::Type) = SparseSymmetricMatrix{T}(Dict{(Int,Int),T}(), 0)
SparseSymmetricMatrix() = SparseSymmetricMatrix(Float64)

function setindex!{T<:Number}(m::SparseSymmetricMatrix{T}, v::T, i::Int, j::Int)
    if i <= j
        if j > m.size
            m.size = j
        end
        m.entries[(i, j)] = v
    else
        if i > m.size
            m.size = i
        end
        m.entries[(j, i)] = v
    end
end

function getindex{T<:Number}(m::SparseSymmetricMatrix{T}, i::Int, j::Int)
    if i <= j
        get(m.entries, (i, j), zero(T))
    else
        get(m.entries, (j, i), zero(T))
    end
end

start(m::SparseSymmetricMatrix) = start(m.entries)
done(m::SparseSymmetricMatrix, state) = done(m.entries, state)
next(m::SparseSymmetricMatrix, state) = next(m.entries, state)



type SparseSymmetricBlockMatrix{T<:Number}
    blocks::Dict{Int, SparseSymmetricMatrix{T}}
end

SparseSymmetricBlockMatrix(T::Type) = SparseSymmetricBlockMatrix(Dict{Int, SparseSymmetricMatrix{T}}())
SparseSymmetricBlockMatrix() = SparseSymmetricBlockMatrix(Float64)

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, v::T, bi::Int, i::Int, j::Int)
    if has(bm.blocks, bi)
        bm.blocks[bi][i, j] = v
    else
        m = SparseSymmetricMatrix(T)
        m[i, j] = v
        bm.blocks[bi] = m
    end
end

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, m::SparseSymmetricMatrix{T}, bi::Int)
    bm.blocks[bi] = m
end

function getindex{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, bi::Int, i::Int, j::Int)
    get(bm.blocks, bi, SparseSymmetricMatrix(T))[i, j]
end

start(m::SparseSymmetricBlockMatrix) = start(m.blocks)
done(m::SparseSymmetricBlockMatrix, state) = done(m.blocks, state)
next(m::SparseSymmetricBlockMatrix, state) = next(m.blocks, state)






type SparseSDP{T<:Number}
    c::SparseSymmetricBlockMatrix{T}
    as::Dict{Int, SparseSymmetricBlockMatrix{T}}
    bs::Dict{Int, T}
end

SparseSDP(T::Type) = SparseSDP(SparseSymmetricBlockMatrix(T), 
                                 Dict{Int,SparseSymmetricBlockMatrix{T}}(),
                                 Dict{Int,T}())

SparseSDP() = SparseSDP(Float64)

setc!{T<:Number}(sdp::SparseSDP{T}, c::SparseSymmetricBlockMatrix{T}) = sdp.c = c

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
        a = SparseSymmetricBlockMatrix(T)
        a[bi, i, j] = v
        sdp.as[ri] = a
    end
end

geta(sdp::SparseSDP, ri::Int) = sdp.tcs[ri]

geta(sdp::SparseSDP, ri::Int, bi::Int) = sdp

setb!{T<:Number}(sdp::SparseSDP{T}, ri::Int, v::T) = sdp.bs[ri] = v

getb(sdp::SparseSDP, ri::Int) = sdp.bs[ri]

function getbvec{T<:Number}(sdp::SparseSDP{T})
    l = zeros(T, m(sdp))
    for (ri, v) in sdp.bs
        l[ri] = v
    end
    l
end
    






function addblocksizes(bsizes::Dict{Int,Int}, m::SparseSymmetricBlockMatrix)
    for (bi, block) in m.blocks
        if get(bsizes, bi, 0) <= block.size
            bsizes[bi] = block.size
        end
    end
    bsizes
end

blocksizes(m::SparseSymmetricBlockMatrix) = addblocksizes(Dict{Int,Int}(), m)

function blocksizes(sdp::SparseSDP)
    bsizes = blocksizes(sdp.c)
    for bm in values(sdp.as)
        addblocksizes(bsizes, bm)
    end
    bsizesvector = zeros(Int, max(0, keys(bsizes)...))
    for (bi, bsize) in bsizes
        bsizesvector[bi] = bsize
    end
    bsizesvector
end

function m(sdp::SparseSDP)
    v = 0
    for ri in keys(sdp.as)
        if ri > v
            v = ri
        end
    end
    for ri in keys(sdp.bs)
        if ri > v
            v = ri
        end
    end
    v
end

function writesdpasparse(io::IO, sdp::SparseSDP)
    bsizes = blocksizes(sdp)    

    println(io, m(sdp))
    println(io, length(bsizes))
    
    for b in bsizes[1:end-1]
        print(io, b)
        print(io, ' ')
    end
    if length(bsizes) > 0
        println(io, bsizes[end])
    end
    
    bvec = getbvec(sdp)
    for r in bvec[1:end-1]
        print(io, r)
        print(io, ' ')
    end
    if length(bvec) > 0
        println(io, bvec[end])
    end
    
    for (bi, block) in sdp.c
        for ((i, j), v) in block
            print(io, '0')
            print(io, ' ')
            print(io, bi)
            print(io, ' ')
            print(io, i)
            print(io, ' ')
            print(io, j)
            print(io, ' ')
            print(io, v)
            println(io)
        end
    end
    
    for (ri, bm) in sdp.as
        for (bi, block) in bm
            for ((i, j), v) in block
                print(io, ri)
                print(io, ' ')
                print(io, bi)
                print(io, ' ')
                print(io, i)
                print(io, ' ')
                print(io, j)
                print(io, ' ')
                print(io, v)
                println(io)
            end
        end
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
    println(datafname)
    println(sdp.as)
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