using SemidefiniteProgramming

sdp = SparseSDP()

setc!(sdp, 1, 1, 1, 1.0)
setc!(sdp, 1, 2, 2, 2.0)
setc!(sdp, 2, 1, 1, 3.0)
setc!(sdp, 2, 2, 2, 4.0)

setb!(sdp, 1, 10.0)
seta!(sdp, 1, 1, 1, 1, 1.0)
seta!(sdp, 1, 1, 2, 2, 1.0)

setb!(sdp, 2, 20.0)
seta!(sdp, 2, 2, 1, 1, 5.0)
seta!(sdp, 2, 2, 1, 2, 2.0)
seta!(sdp, 2, 2, 2, 2, 6.0)

println(solve(sdp, CSDP()))
