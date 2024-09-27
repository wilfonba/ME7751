# Compiling
This code is compiled by running `./make.sh` in the directory this README.md is
in.

# Running
The problem parameters are input using the namelist file `main.inp`.
The input options are:
- gamma0: The coefficient for the zero degree term in the polynomial
          defining gamma
- gamma1: The coefficient for the first degree term in the polynomial
          defining gamma
- L: The length of the domain
- N: The number of cells to use
- Q0: The coefficient for the zero degree term in the polynomial
      defining Q
- Q1: The coefficient for the first degree term in the polynomial
      defining Q
- U: The convective velocity
- tol: Convergence tolerance for nonlinear problems
- bench: Toggles benchmarking mode. When true, the code is ran 10^6
         times and the average runtime is output.

# Output
Before performing the finite difference solve, the code writes the
following to the terminal:
```
 Solving the 1D Steady Transport Equation with:
 gamma0:   0.10000000000000001
 gamma1:   0.10000000000000001
 U:    0.0000000000000000
 Q0:    0.0000000000000000
 Q1:    0.0000000000000000
 N:           20
```
This serves as sanity check for the user inputs. After execution,
the code writes the following to the terminal:
```
 Execution time:    4.99999987E-05  s
 Iterations:           12
```
This gives information about runtime and, when the problem is nonlinear,
the number of iterations required to reach convergence. The solution is
written to the file `output.csv`. The first column of this file contains
the x coordinates of the finite difference points and the second column
contains the solution value.
