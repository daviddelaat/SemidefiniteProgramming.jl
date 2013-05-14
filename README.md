## Semidefinite Programming

This package provides a Julia interface to model semidefinite programming problems and solve them using solvers such as SDPA and CSDP. The library is in development, so many features are incomplete and changes (also to the API) are expected.

# Introduction

Consider a semidefinite program of the form

    max tr(C X) subject to X is positive semidefinite
                           tr(A_i X) = b_i for i = 1, ...,m

Here `C`, `X`, and `A_1`, `...`, `A_m` are symmetric block matrices (all assumed to have identical size and block structure), and `b_1`, `...`, `b_m` are scalars. 

This problem can be modelled by constructing a sparse semidefinite program with

```julia
    using SemidefiniteProgramming
    sdp = SparseSDP(maximize=true) 
```

and then setting the nonzero scalars and the nonzero entries of the matrices. The most basic way to do this is as follows: For the scalars `b_i` use

```julia
    setrhs!(sdp, i, value)
```

For the entries of the objective matrix `C` use

```julia
    setobj!(sdp, blockindex, rowindex, columnindex, value)
```

For the constraint matrices `A_i` use

```julia
    setcon!(sdp, i, blockindex, rowindex, columnindex, value)
```

Then we solve the program with 

```julia
    solve(sdp, solver)
```

Notice that the number of constraints, the number of blocks, and the blocksizes do not need to be specified; they will be determined automatically  based on the entries that are set. Of course all the matrices involded are assumed to have identical block structure. The indices of the contraints, blocks, and matrices do not need be integers; all Julia object can be used here. When storing a SparseSDP in for instance the SDPA-sparse format the indices will be converted to integers.

# Example

Consider the program above with `b1 = 10`, `b2 = 20`, 

```
C = [1 0 0 0;
     0 2 0 0;
     0 0 3 0;
     0 0 0 4],
```

```
A1 = [1 0 0 0;
      0 1 0 0;
      0 0 0 0;
      0 0 0 0],
```

and

    A2 = [0 0 0 0;
          0 1 0 0;
          0 0 5 2;
          0 0 2 6]

To solve this we use

```julia
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
setcon!(sdp, 2, 1, 1, 1, 1.0)
setcon!(sdp, 2, 2, 1, 1, 1.0)
setcon!(sdp, 2, 3, 1, 1, 5.0)
setcon!(sdp, 2, 3, 1, 2, 2.0)
setcon!(sdp, 2, 3, 2, 2, 6.0)

println(solve(sdp, CSDP()))
```