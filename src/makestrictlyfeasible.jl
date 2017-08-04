function makeprimalstrictlyfeasible!{T<:Number}(sdp::SparseSDP{T}, epsilon::Number)
    for (ci, v) in rhs(sdp)
        v = zero(T)
        blockmatrix = cons(sdp)[ci]
        for (_, block) in blockmatrix
            for (indices, entry) in block
                i, j = indices
                if i == j
                    v += entry
                end
            end
        end
        setrhs!(sdp, ci, rhs(sdp, ci) + epsilon*v)
    end
end

function makedualstrictlyfeasible!{T<:Number}(sdp::SparseSDP{T}, epsilon::Number)
    for (_, block) in obj(sdp)
        for (indices, value) in block
            i, j = indices
            if i == j
                block[i, j] = block[i,j] - epsilon
            end
        end
    end 
end
