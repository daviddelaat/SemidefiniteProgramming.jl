using SemidefiniteProgramming

sdp = SparseSDP()

setobj!(sdp, 1, 1, 1, 1.0)
setobj!(sdp, 1, 2, 2, 2.0)
setobj!(sdp, "z", :x, :x, 3.0)
setobj!(sdp, "z", 2, 2, 4.0)

setrhs!(sdp, 0.0, 10.0)
setcon!(sdp, 0.0, 1, 1, 1, 1.0)
setcon!(sdp, 0.0, 1, 2, 2, 1.0)

setrhs!(sdp, 2, 20.0)
setcon!(sdp, 2, "z", :x, :x, 5.0)
setcon!(sdp, 2, "z", :x, 2, 2.0)
setcon!(sdp, 2, "z", 2, 2, 6.0)

println(blocksizes(sdp))

println(solve(sdp, CSDP()))
