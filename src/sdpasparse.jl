function writesdpasparse(io::IO, sdp::SparseSDP)
    if !normalized(sdp)
        sdp, = normalize(sdp)
    end

    if !ismaximizationproblem(sdp)
        sdp.obj = -1.0 * sdp.obj
    end
    
    bsizes = blocksizes(sdp)

    println(io, ncons(sdp))
    println(io, length(bsizes))
    for i = 1:length(bsizes)
        print(io, bsizes[i])
        print(io, ' ')
    end
    println(io)
    
    rhs = rhsdense(sdp)
    for r in rhs[1:end-1]
        print(io, r)
        print(io, ' ')
    end
    if length(rhs) > 0
        println(io, rhs[end])
    end
    
    for (bi, block) in obj(sdp)
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
    
    for (ci, bm) in cons(sdp)
        for (bi, block) in bm
            for ((i, j), v) in block
                print(io, ci)
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
    sdp = NormalizedSparseSDP()
    
    s = 1
    for l in readlines(io)
        if !beginswith(l, '#')
            if s == 4
                t = split(l) 
                for i = 1:length(t)
                     setrhs!(sdp, i, parsefloat(t[i]))
                end
            elseif s > 4
                t = split(l)
                ci = parseint(t[1])
                bi = parseint(t[2])
                i = parseint(t[3])
                j = parseint(t[4])
                v = parsefloat(t[5])
                if ci == 0
                    setobj!(sdp, bi, i, j, v)
                else
                    setcon!(sdp, ci, bi, i, j, v)
                end
            end 
            s += 1
        end
    end
    
    sdp
end
