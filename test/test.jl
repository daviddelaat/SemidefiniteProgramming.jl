using SemidefiniteProgramming

sdp = SparseSDP()

setc!(sdp, 1, 1, 1, 1.0)
setc!(sdp, 1, 2, 2, 2.0)
setc!(sdp, "z", :x, :x, 3.0)
setc!(sdp, "z", 2, 2, 4.0)

setb!(sdp, 0.0, 10.0)
seta!(sdp, 0.0, 1, 1, 1, 1.0)
seta!(sdp, 0.0, 1, 2, 2, 1.0)

setb!(sdp, 2, 20.0)
seta!(sdp, 2, "z", :x, :x, 5.0)
seta!(sdp, 2, "z", :x, 2, 2.0)
seta!(sdp, 2, "z", 2, 2, 6.0)

println(blocksizes(sdp))

println(solve(sdp, CSDP()))
