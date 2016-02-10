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



function solve{T<:Real}(sdp::SparseSDP{T}, solver::SDPAGEN; io::IO=STDOUT, outputfname="output.temp")
    datafname, dataio = mktemp()
    if verbose(solver)
        println(io, "Data filename: $datafname")
    end
    writesdpasparse(dataio, sdp)
    flush(dataio)
    primalobjective = NaN
    dualobjective = NaN
   
    for l in eachline(`$(executable(solver)) -ds $datafname -o $outputfname`)
        if verbose(solver)
            print(io, l)
        end
        if startswith(l, "objValPrimal = ")
            f = float(strip(split(l, " = ")[2]))
            primalobjective = ismaximizationproblem(sdp) ? f : -f
        end
        if startswith(l, "objValDual ")
            f = float(strip(split(l, " = ")[2]))
            dualobjective = ismaximizationproblem(sdp) ? f : -f
         end
    end
    close(dataio)
    SparseSDPSolution(primalobjective, dualobjective) 
end

function solve(sdp::SparseSDP, solver::CSDP; io::IO=STDOUT, ioverbose::IO=STDOUT, outputfname=nothing)
    sdp, cm, bm, ems = normalize(sdp)

    datafname, dataio = mktemp()
    if outputfname == nothing
        outputfname, outputio = mktemp()
    else
        outputio = open(outputfname, "w")
    end
    if verbose(solver)
        println(io, "Data filename: $datafname")
        println(io, "Output filename: $outputfname")
    end
    writesdpasparse(dataio, sdp)
    close(dataio)
    
    primalobjective = NaN
    dualobjective = NaN
    for l in eachline(`$(executable(solver)) $datafname $outputfname`)
        print(io, l)
        flush(io)
        if verbose(solver)
            print(ioverbose, l)
        end
        if startswith(l, "Primal objective value: ")
            f = float(strip(split(l, ": ")[2]))
            primalobjective = ismaximizationproblem(sdp) ? f : -f
        end
        if startswith(l, "Dual objective value: ")
            f = float(strip(split(l, ": ")[2]))
            dualobjective = ismaximizationproblem(sdp) ? f : -f
        end
    end
    sol = SparseSDPSolution(primalobjective, dualobjective)
    readcsdpoutput!(outputio, sol, cm, bm, ems)
    close(outputio)
    sol   
end
