type IndexMap
    m::Dict
    c::Int
end

IndexMap() = IndexMap(Dict(), 0)

function newindex!(oldindex, i::IndexMap)
    if haskey(i.m, oldindex)
        return i.m[oldindex]
    else
        i.c += 1
        i.m[oldindex] = i.c
        return i.c        
    end
end

function normalize(sdp::SparseSDP)
    newsdp = SparseSDP(maximize=ismaximizationproblem(sdp))
    
    cm = IndexMap()
    bm = IndexMap()
    em = Dict()
    
    for (ri, value) in rhs(sdp)
        newri = newindex!(ri, cm)
        setrhs!(newsdp, newri, value)
    end
    
    for (bi, matrix) in blocks(obj(sdp))
        newbi = newindex!(bi, bm)
        if !haskey(em, newbi)
            em[newbi] = IndexMap()
        end
        for ((i, j), v) in entries(matrix)
            newi = newindex!(i, em[newbi])
            newj = newindex!(j, em[newbi])
            setobj!(newsdp, newbi, newi, newj, v)
        end    
    end
       
    for (ri, blockmatrix) in cons(sdp)
        newri = newindex!(ri, cm)
        for (bi, matrix) in blocks(blockmatrix)
            newbi = newindex!(bi, bm)
            if !haskey(em, newbi)
                em[newbi] = IndexMap()
            end
            for ((i, j), v) in entries(matrix)
                newi = newindex!(i, em[newbi])
                newj = newindex!(j, em[newbi])
                setcon!(newsdp, newri, newbi, newi, newj, v)
            end    
        end
    end

    newsdp
end
