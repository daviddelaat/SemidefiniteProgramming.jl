abstract SDPSolver

verbose(solver::SDPSolver) = solver.verbose


abstract SDPAGEN <: SDPSolver

executable(solver::SDPAGEN) = solver.executable


immutable SDPA <: SDPAGEN
    speed::Int
    executable::ASCIIString
    verbose::Bool
end

SDPA(; speed::Integer=0, 
       executable::ASCIIString="sdpa", 
       verbose::Bool=false) = SDPA(speed, executable, verbose)


immutable SDPAQD <: SDPAGEN
    speed::Int
    executable::ASCIIString
    verbose::Bool
end

SDPAQD(; speed::Integer=0, 
         executable::ASCIIString="sdpa_qd", 
         verbose::Bool=false) = SDPAQD(speed, executable, verbose)


immutable SDPAGMP <: SDPAGEN
    speed::Int
    executable::ASCIIString
    verbose::Bool
end

SDPAGMP(; speed::Integer=0, 
          executable::ASCIIString="sdpa_gmp", 
          verbose::Bool=false) = SDPAGMP(speed, executable, verbose)


immutable CSDP <: SDPSolver
    executable::ASCIIString
    verbose::Bool
end

CSDP(; executable="csdp", verbose=false) = CSDP(executable, verbose)

executable(solver::CSDP) = solver.executable



function solve(sdp::SparseSDP, solver::SDPAGEN)
    datafname, dataio = mktemp()
    if verbose(solver)
        println("Data filename: $datafname")
    end
    writesdpasparse(dataio, sdp)
    flush(dataio)
    for l in eachline(`$(solver.executable) -ds $datafname -o /dev/null`)
        if verbose(solver)
            println(l)
        end
        if beginswith(l, "objValPrimal = ")
            close(dataio)
            f = float(split(l, " = ")[2])
            return ismaximizationproblem(sdp) ? f : -f
        end      
    end
    close(dataio)
    nothing
end

function solve(sdp::SparseSDP, solver::CSDP)
    datafname, dataio = mktemp()
    if verbose(solver)
        println("Data filename: $datafname")
    end
    writesdpasparse(dataio, sdp)
    flush(dataio)    
    for l in eachline(`$(solver.executable) $datafname /dev/null`)
        if verbose(solver)
            println(l)
        end
        if beginswith(l, "Primal objective value: ")
            close(dataio)
            f = float(split(l, ": ")[2])
            return ismaximizationproblem(sdp) ? f : -f
        end
    end
    close(dataio)
    nothing
end