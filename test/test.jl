using SemidefiniteProgramming

sdp = SparseSDP()

setobj!(sdp, 1, 1, 2, 1.0)

setcon!(sdp, 1, 1, 1, 1, 1.0)
setrhs!(sdp, 1, 1.0)

setcon!(sdp, 2, 1, 2, 2, 1.0)
setrhs!(sdp, 2, 1.0)

sol = solve(sdp, CSDP())

println(primalmatrix(sol))
