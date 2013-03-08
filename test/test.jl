using SemidefiniteProgramming

sdp = SparseSDP([10., 20.], 
                [Rel(0, 1, 1, 1, 1.),
                 Rel(0, 1, 2, 2, 2.),
                 Rel(0, 2, 1, 1, 3.),
                 Rel(0, 2, 2, 2, 4.),
                 Rel(1, 1, 1, 1, 1.),
                 Rel(1, 1, 2, 2, 1.),
                 Rel(2, 1, 2, 2, 1.),
                 Rel(2, 2, 1, 1, 5.),
                 Rel(2, 2, 1, 2, 2.),
                 Rel(2, 2, 2, 2, 6.)])

println(solve(sdp, CSDP()))



