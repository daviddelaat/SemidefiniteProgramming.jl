using SemidefiniteProgramming

sdp = SparseSDP(maximize=true)

setobj!(sdp, 1, 1, 1, 1.0)
setobj!(sdp, 2, 1, 1, 2.0)
setobj!(sdp, 3, 1, 1, 3.0)
setobj!(sdp, 3, 2, 2, 4.0)

setrhs!(sdp, 1, 10.0)
setcon!(sdp, 1, 1, 1, 1, 1.0)
setcon!(sdp, 1, 2, 1, 1, 1.0)

setrhs!(sdp, 2, 20.0)
setcon!(sdp, 2, 2, 1, 1, 1.0)
setcon!(sdp, 2, 3, 1, 1, 5.0)
setcon!(sdp, 2, 3, 1, 2, 2.0)
setcon!(sdp, 2, 3, 2, 1, 2.0)
setcon!(sdp, 2, 3, 2, 2, 6.0)

println(obj(solve(sdp, CSDP())))

