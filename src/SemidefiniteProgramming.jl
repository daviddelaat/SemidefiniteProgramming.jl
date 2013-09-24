module SemidefiniteProgramming

using Cartesian

import 
    Base.print, 
    Base.show,
    Base.getindex,
    Base.setindex!,
    Base.size,
    Base.start,
    Base.next,
    Base.done,
    Base.copy,
    Base.(*),
    Base.delete!
    
export 
    SparseSDP,
    ismaximizationproblem,
    setmaximizationproblem!,
    setobj!,
    obj,
    setrhs!,
    rhs,
    setcon!,
    cons,
    ncons,
    blocksizes,
    
    writesdpasparse,
    readsdpasparse,
    
    SDPSolver,
    SDPA,
    SDPAQD,
    SDPAGMP,
    CSDP,
    solve,

    addsoscon!

include("sparsesymmetricmatrix.jl")
include("sparsesymmetricblockmatrix.jl")
include("sparsesdp.jl")
include("normalization.jl")
include("sdpasparse.jl")
include("solvers.jl")
include("sos.jl")

end
