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

    for ri in keys(rhs(sdp))
        if isa(ri, Int)
	        add!(integerindices, ri)
        else
	        add!(strangeindices, ri)
        end
    end

    for ri in keys(cons(sdp))
        if isa(ri, Int)
	        add!(integerindices, ri)
        else
            add!(strangeindices, ri)
        end
    end

    for strangeindex in strangeindices
        integerindex = nextavailable(integerindices)
	    used = false
	    if haskey(rhs(sdp), strangeindex)
	        v = rhs(sdp, strangeindex)
	        delete!(rhs(sdp), strangeindex)
	        setrhs!(sdp, integerindex, v)
	        used = true
	    end
	    if haskey(cons(sdp), strangeindex)
    	    m = getcon(sdp, strangeindex)
	        delete!(cons(sdp), strangeindex)
	        cons(sdp)[integerindex] = m
	        used = true
	    end
	    if used
	        add!(integerindices, integerindex)
	    end
    end 
end

function normalize_matrix_indices!(sdp::SparseSDP)
    for blockindex in keys(blocksizes(sdp))
        ms = Array(SparseSymmetricMatrix, length(cons(sdp)) + 1)
        ms[1] = obj(sdp)[blockindex]
        ms[2:end] = [x[blockindex] for x in values(cons(sdp))]
        
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
     
    for blockindex in keys(blocks(obj(sdp)))
        add!(isa(blockindex, Int) ? integerindices : strangeindices, blockindex)
    end
     
    for blockmatrix in values(cons(sdp))
       for blockindex in keys(blocks(blockmatrix))
           add!(isa(blockindex, Int) ? integerindices :strangeindices, blockindex)
       end
    end
    
    for strangeindex in strangeindices
        if haskey(blocks(obj(sdp)), strangeindex)
            v = blocks(obj(sdp))[strangeindex]
            delete!(blocks(obj(sdp)), strangeindex)
            blocks(obj(sdp))[nextavailable(integerindices)] = v
        end
        for blockmatrix in values(cons(sdp))
            if haskey(blocks(blockmatrix), strangeindex)
                v = blocks(blockmatrix)[strangeindex]
                delete!(blocks(blockmatrix), strangeindex)
                blocks(blockmatrix)[nextavailable(integerindices)] = v
            end
        end
    end
end

function normalize_indices!(sdp::SparseSDP)
    normalize_block_indices!(sdp)
    normalize_constraint_indices!(sdp)
    normalize_matrix_indices!(sdp)
end
