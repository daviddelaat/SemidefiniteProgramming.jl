abstract SDPSolver

verbose(solver::SDPSolver) = solver.verbose

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
    verbose::Bool
end

CSDP(; executable="csdp", verbose=false) = CSDP(executable, verbose)

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
    sdp = copy(sdp)
    if !ismaximizationproblem(sdp)
        obj(sdp) = -obj(sdp)
    end
    normalize_indices!(sdp)
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