using SemidefiniteProgramming

sdp = SparseSDP()

sdp[0][1][1,1] = 1.0
sdp[0,1,2,2] = 2.0
sdp[0,2,1,1] = 3.0
sdp[0,2,2,2] = 4.0

sdp[1,0] = 10.0
sdp[1,1,1,1] = 1.0
sdp[1,1,2,2] = 1.0

sdp[2,0] = 20.0
sdp[2,2,1,1] = 5.0
sdp[2,2,1,2] = 2.0
sdp[2,2,2,2] = 6.0

println(sdp)

println(solve(sdp, CSDP()))
