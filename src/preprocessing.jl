function get_equal_constraints(sdp::SparseSDP)
    pairs = Set()
    indices = collect(keys(cons(sdp)))
    for i = 1:length(indices), ip = i+1:length(indices)
        ci = indices[i]
        cip = indices[ip]
        bm = cons(sdp)[ci]
        bmp = cons(sdp)[cip]
        if ci != cip && isapprox(bm, bmp) && isapprox(get(rhs(sdp), ci, 0.0), get(rhs(sdp), cip, 0.0))
            push!(pairs, (ci, cip))
        end
    end
    pairs
end

function remove_equal_constraints!(sdp::SparseSDP)
    eqs = get_equal_constraints(sdp)
    for (ci, _) in eqs
         delete!(cons(sdp), ci)
         delete!(rhs(sdp), ci)
     end
end

function find_removable_row(A::Matrix)
    m = size(A, 1)
    n = size(A, 2)
    r = rank(A)
    for i = 1:m
        B = copy(A)
        B[i, :] = 0.0
        if rank(B) == r
            return i
        end
    end
end

function remove_linear_dependence!(sdp::SparseSDP; debug::Bool=false, dorref=false)
    constraintindices = union(keys(rhs(sdp)), keys(cons(sdp)))
    m = length(constraintindices)

    elindices = Dict()
    n = 1
    debug && println("building dictionaries")
    for bm in values(cons(sdp))
        for (bi, b) in blocks(bm)
            for ((ri, ci), _) in b
                if !haskey(elindices, (bi,ri,ci))
                    elindices[(bi, ri, ci)] = n 
                    n += 1
                end
            end
        end
    end
    
    debug && println("building matrix")
    A = zeros(n, m)
    for j = 1:m
        A[n, j] = get(rhs(sdp), constraintindices[j], 0.0)
        for (bi, b) in cons(sdp)[constraintindices[j]]
            for ((ri, ci), v) in b
                A[elindices[(bi, ri, ci)], j] = v
            end
        end
    end
    
#        v = []
#        for j = 1:m
#            w = vcat(v, [j])
#            if rank(A[:, w]) == length(w)
#                v = w
#            else
#                ci = constraintindices[j]
#                delete!(rhs(sdp), ci)
#                delete!(cons(sdp), ci)
#            end
#        end
#    println(size(A))
#    U,S,V = svd(A)
#    println(S)
#    debug && println("removing constraints")
#    if dorref
#        debug && println("rref")    
#        @time B = rref(A)
#        i = 1
#        counter=0
#        for j = 1:m
#            if isapprox(B[i, j], 0.0, atol=10e-8)
#                counter+=1
#                ci = constraintindices[j]
#                delete!(rhs(sdp), ci)
#                delete!(cons(sdp), ci)
#            else
#                i += 1
#            end
#        end
#        println(counter)
#    else
    println(size(A))
    _, r, p = qr(A, Val{true})
    
    for i = size(r, 1):-1:1
        println(r[i,i])
        if abs(r[i, i]) > 1e-14
            break
        end
        ci = constraintindices[p[i]]
        delete!(rhs(sdp), ci)
        delete!(cons(sdp), ci)
    end
#        i = 1
#        counter = 0
#        for j = 1:size(r, 2)
#            println(i, " " ,  j)
#            if abs(r[i,j]) >= 10e-10
#                i += 1
#            else
#                println(r[i,j])
#                counter += 1
#                ci = constraintindices[p[j]]
#                delete!(rhs(sdp), ci)
#                delete!(cons(sdp), ci)
#            end
#            if i > size(r, 1)
#                break
#            end
#        end
#        println(counter)
end

