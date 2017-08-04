mutable struct IndexMap
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

function normalize{T}(sdp::SparseSDP{T})
    newsdp = SparseSDP(T, maximize=ismaximizationproblem(sdp), normalized=true)
    
    cm = IndexMap()
    bm = IndexMap()
    ems = Dict{Any,IndexMap}()
    
    for (ri, value) in rhs(sdp)
        newri = newindex!(ri, cm)
        setrhs!(newsdp, newri, value)
    end
    
    for (bi, matrix) in blocks(obj(sdp))
        newbi = newindex!(bi, bm)
        if !haskey(ems, newbi)
            ems[newbi] = IndexMap()
        end
        for ((i, j), v) in entries(matrix)
            newi = newindex!(i, ems[newbi])
            newj = newindex!(j, ems[newbi])
            setobj!(newsdp, newbi, newi, newj, v)
        end    
    end
       
    for (ri, blockmatrix) in cons(sdp)
        newri = newindex!(ri, cm)
        for (bi, matrix) in blocks(blockmatrix)
            newbi = newindex!(bi, bm)
            if !haskey(ems, newbi)
                ems[newbi] = IndexMap()
            end
            for ((i, j), v) in entries(matrix)
                newi = newindex!(i, ems[newbi])
                newj = newindex!(j, ems[newbi])
                setcon!(newsdp, newri, newbi, newi, newj, v)
            end
        end
    end
    nems = Dict()
    for (k,v) in ems
        nems[k] = ems[k].m
    end
    newsdp, cm.m, bm.m, nems
end
