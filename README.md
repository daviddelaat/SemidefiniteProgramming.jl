# Semidefinite Programming

This package provides a Julia interface for low-level modeling of semidefinite programming problems and for solving semidefinite programs with solvers such as SDPA and CSDP.

Maintenance status: Currently no new futures are being developed for this package. Bugs will be fixed and this package will be kept up to date with new Julia releases.

# Introduction

Consider a semidefinite program of the form

    max tr(C X) subject to X is positive semidefinite
                           tr(A_i X) = b_i for i = 1, ...,m

Here `C`, `X`, and `A_1`, `...`, `A_m` are symmetric block matrices (all assumed to have identical size and block structure), and `b_1`, `...`, `b_m` are scalars. 

This problem can be modeled by constructing a sparse semidefinite program with

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
sol = solve(sdp, solver)
```
and print the (primal) objective value:
```julia
println(obj(sol))
```

Notice that the number of constraints, the number of blocks, and the blocksizes do not need to be specified; they will be determined automatically based on the entries you have set. Of course all the matrices involded are assumed to have identical block structure. The indices of the contraints, blocks, and matrices do not need be integers; you can use any Julia object here. When storing a SparseSDP in for instance the SDPA-sparse format the indices will be converted to integers automatically. 

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
setcon!(sdp, 2, 2, 1, 1, 1.0)
setcon!(sdp, 2, 3, 1, 1, 5.0)
setcon!(sdp, 2, 3, 1, 2, 2.0)
setcon!(sdp, 2, 3, 2, 1, 2.0)
setcon!(sdp, 2, 3, 2, 2, 6.0)

println(obj(solve(sdp, CSDP())))
```

# Solvers

To use a solver construct an immutable solver object with `CSDP()`, `SDPA()`, etc, and supply it as the second argument to the `solve` function. The solver objects support the optional named arguments 
 - `verbose` (print solver output to stdout)
 - `executable` (path to the solver executable)

## CSDP

To use the CSDP solver you need to install [CSDP](https://projects.coin-or.org/Csdp) and make sure that the CSDP binary is in your path. On Debian/Ubuntu you can do this by installing the `coinor-csdp` package. On Fedora it is the `csdp` package. 

## SDPA

To use one of the SDPA solvers install [SDPA](http://sdpa.sourceforge.net/) and make sure the executable is in your path. On Debian/Ubuntu you can do this by installing the package `sdpa` (this package only contains the standard SDPA solver). Use SDPA for the standard SDPA solver and SDPAQD or SDPAGMP for the high precision solvers.

# SparseSDPSolution objects

Having solved a semidefinite program with
```julia
sol = solve(sdp, CSDP())
```
you can extract the primal and dual objective values with `obj(sol)` and `dualobj(sol)`. To extract the values of the optimal primal variables (the matrix `X` in the notation above) use
```julia
primalmatrix(sol)[blockindex][rowindex, columnindex]
```
Variable extraction is currently only supported with the CSDP solver.
