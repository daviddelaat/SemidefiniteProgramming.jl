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
    writesdpasparse,
    readsdpasparse,
    SDPSolver,
    SDPA,
    SDPAQD,
    SDPAGMP,
    CSDP,
    solve

   
    

    
type SparseSymmetricMatrix{T<:Number}
    entries::Dict{(Any, Any), T}
    indices::Set
end

SparseSymmetricMatrix(T::Type) = SparseSymmetricMatrix{T}(Dict{(Any,Any),T}(), Set())
SparseSymmetricMatrix() = SparseSymmetricMatrix(Float64)

entries(m::SparseSymmetricMatrix) = m.entries
indices(m::SparseSymmetricMatrix) = m.indices

function setindex!{T<:Number}(m::SparseSymmetricMatrix{T}, v::T, i, j)
    if method_exists(isless, (typeof(i), typeof(j)))
        if i <= j
            entries(m)[(i, j)] = v
        else
            entries(m)[(j, i)] = v
        end
    elseif method_exists(hash, (typeof(i),)) && method_exists(hash, (typeof(j),)) 
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

function getindex{T<:Number}(m::SparseSymmetricMatrix{T}, i, j)
    if method_exists(isless, (typeof(i), typeof(j)))
        if i <= j
            get(m.entries, (i, j), zero(T))
        else
            get(m.entries, (j, i), zero(T))
        end
    elseif method_exists(hash, (typeof(i),)) && method_exists(hash, (typeof(j),)) 
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

function smallestavailableindex(ms::Vector{SparseSymmetricMatrix})
    s = 1
    for m in ms
        for (i, j) in keys(m.entries)
            if is(i, Int) && i > s
                s = i
            end
            if is(j, Int) && j > s
                s = j
            end
        end
    end
    s
end




type SparseSymmetricBlockMatrix{T<:Number}
    blocks::Dict{Any,SparseSymmetricMatrix{T}}
end

SparseSymmetricBlockMatrix(T::Type) = SparseSymmetricBlockMatrix(Dict{Any,SparseSymmetricMatrix{T}}())
SparseSymmetricBlockMatrix() = SparseSymmetricBlockMatrix(Float64)

function setindex!{T<:Number}(bm::SparseSymmetricBlockMatrix{T}, v::T, bi, i, j)
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

function getindex(bm::SparseSymmetricBlockMatrix, bi, i, j)
    get(bm.blocks, bi, SparseSymmetricMatrix(T))[i, j]
end

getindex{T}(bm::SparseSymmetricBlockMatrix{T}, bi) = get(bm.blocks, bi, SparseSymmetricMatrix(T))

start(m::SparseSymmetricBlockMatrix) = start(m.blocks)
done(m::SparseSymmetricBlockMatrix, state) = done(m.blocks, state)
next(m::SparseSymmetricBlockMatrix, state) = next(m.blocks, state)






type SparseSDP{T<:Number}
    c::SparseSymmetricBlockMatrix{T}
    as::Dict{Any,SparseSymmetricBlockMatrix{T}}
    bs::Dict{Any,T}
end

SparseSDP(T::Type) = SparseSDP(SparseSymmetricBlockMatrix(T), 
                               Dict{Any,SparseSymmetricBlockMatrix{T}}(),
                               Dict{Any,T}())

SparseSDP() = SparseSDP(Float64)

setc!{T<:Number}(sdp::SparseSDP{T}, c::SparseSymmetricBlockMatrix{T}) = sdp.c = c

function setc!{T<:Number}(sdp::SparseSDP{T}, bi, m::SparseSymmetricMatrix{T})
    sdp.c[bi] = m
end

function setc!{T<:Number}(sdp::SparseSDP{T}, bi, i, j, v::T)
    sdp.c[bi, i, j] = v
end

getc(sdp::SparseSDP) = sdp.c

function seta!{T<:Number}(sdp::SparseSDP{T}, ri, bi, i, j, v::T)
    if has(sdp.as, ri)
        sdp.as[ri][bi, i, j] = v
    else
        a = SparseSymmetricBlockMatrix(T)
        a[bi, i, j] = v
        sdp.as[ri] = a
    end
end

geta(sdp::SparseSDP, ri) = sdp.tcs[ri]

geta(sdp::SparseSDP, ri, bi) = sdp

setb!{T<:Number}(sdp::SparseSDP{T}, ri, v::T) = sdp.bs[ri] = v

getb(sdp::SparseSDP, ri) = sdp.bs[ri]

function getbvec{T<:Number}(sdp::SparseSDP{T})
    l = zeros(T, m(sdp))
    for (ri, v) in sdp.bs
        l[ri] = v
    end
    l
end






function blocksizes(sdp::SparseSDP)
    blockdict = Dict{Any,Set}()
    for (bi, block) in sdp.c
        indices = get(blockdict, bi, Set())
        add_each!(indices, block.indices)
        blockdict[bi] = indices
    end
    for a in values(sdp.as)
        for (bi, block) in a
            indices = get(blockdict, bi, Set())
            add_each!(indices, block.indices)
            blockdict[bi] = indices
        end
    end
    blocksizedict = Dict{Any,Int}()
    for (bi, entries) in blockdict
        blocksizedict[bi] = length(entries)
    end
    blocksizedict
end


 
function nextavailable(s::Set{Int})
    i = 1
    while contains(s, i)
       i += 1
    end
    i
end





function normalize_constraint_indices!(sdp::SparseSDP)
    integerindices = Set{Int}()
    strangeindices = Set()

    for ri in keys(sdp.bs)
        if isa(ri, Int)
	        add!(integerindices, ri)
        else
	        add!(strangeindices, ri)
        end
    end

    for ri in keys(sdp.as)
        if isa(ri, Int)
	        add!(integerindices, ri)
        else
            add!(strangeindices, ri)
        end
    end

    for strangeindex in strangeindices
        integerindex = nextavailable(integerindices)
	    used = false
	    if has(sdp.bs, strangeindex)
	        v = sdp.bs[strangeindex]
	        delete!(sdp.bs, strangeindex)
	        sdp.bs[integerindex] = v
	        used = true
	    end
	    if has(sdp.as, strangeindex)
    	    m = sdp.as[strangeindex]
	        delete!(sdp.as, strangeindex)
	        sdp.as[integerindex] = m
	        used = true
	    end
	    if used
	        add!(integerindices, integerindex)
	    end
    end 
end

function normalize_matrix_indices!(sdp::SparseSDP)
    for blockindex in keys(blocksizes(sdp))
        ms = Array(SparseSymmetricMatrix, length(sdp.as) + 1)
        ms[1] = sdp.c[blockindex]
        ms[2:end] = [x[blockindex] for x in values(sdp.as)]
        
        integerindices = Set{Int}()
        strangeindices = Set()
    
        for m in ms
            for i in indices(m)
                add!(isa(i, Int) ? integerindices : strangeindices, i)
            end
        end
    
        for strangeindex in strangeindices
            integerindex = nextavailable(integerindices)
            add!(integerindices, integerindex)
            for m in ms
                for ((i, j), v) in entries(m)
                    if i == j == strangeindex
                        delete!(m.entries, (i, j))
                        m[integerindex, integerindex] = v
                    elseif i == strangeindex
                        delete!(m.entries, (i, j))
                        m[integerindex, j] = v
                    elseif j == strangeindex
                        delete!(m.entries, (i, j))
                        m[i, integerindex] = v
                    end
                end
            end
        end             
    end
end
    
function normalize_block_indices!(sdp::SparseSDP)   
    strangeindices = Set()
    integerindices = Set{Int}()
     
    for blockindex in keys(sdp.c.blocks)
        add!(isa(blockindex, Int) ? integerindices : strangeindices, blockindex)
    end
     
    for blockmatrix in values(sdp.as)
       for blockindex in keys(blockmatrix.blocks)
           add!(isa(blockindex, Int) ? integerindices :strangeindices, blockindex)
       end
    end
    
    for strangeindex in strangeindices
        if has(sdp.c.blocks, strangeindex)
            v = sdp.c.blocks[strangeindex]
            delete!(sdp.c.blocks, strangeindex)
            sdp.c.blocks[nextavailable(integerindices)] = v
        end
        for blockmatrix in values(sdp.as)
            if has(blockmatrix.blocks, strangeindex)
                v = blockmatrix.blocks[strangeindex]
                delete!(blockmatrix.blocks, strangeindex)
                blockmatrix.blocks[nextavailable(integerindices)] = v
            end
        end
    end
end

function normalize_indices!(sdp::SparseSDP)
    normalize_block_indices!(sdp)
    normalize_constraint_indices!(sdp)
    normalize_matrix_indices!(sdp)
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
    normalize_indices!(sdp)
    bsizes = values(blocksizes(sdp))    

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

function readsdpasparse(io::IO)
    sdp = SparseSDP()
    
    s = 1
    for l in readlines(io)
        if !beginswith(l, '#')
            if s == 4
                t = split(l) 
                for i = 1:length(t)
                     setb!(sdp, i, parsefloat(t[i]))
                end
            elseif s > 4
                t = split(l)
                ri = parseint(t[1])
                bi = parseint(t[2])
                i = parseint(t[3])
                j = parseint(t[4])
                v = parsefloat(t[5])
                if ri == 0
                    setc!(sdp, bi, i, j, v)
                else
                    seta!(sdp, ri, bi, i, j, v)
                end
            end 
            s += 1
        end
    end
    
    sdp
end


function show{T<:Number}(io::IO, sdp::SparseSDP{T})
#    bsizes = blocksizes(sdp)
#    nr = nrels(sdp)
#    
#    print(io, "SparseSDP with blocksizes ")
#    println(io, bsizes)
#    
#    for i in 0:nr
#        println(io, "Relation $i with right hand side $(get(sdp.rhs, i, zero(T))): ")
#        
#    end
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
    normalize_indices!(sdp)
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
