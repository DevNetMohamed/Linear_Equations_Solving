# Linear Equation System Solver in C++

A small but educational C++ console application that parses linear equations written in natural text form, builds the coefficient matrix, and solves the system using **Gaussian elimination** with partial pivoting.

Originally started as a basic project using Cramer's rule (2019), later significantly improved with better parsing, proper matrix handling, and a more stable & efficient solver.

## Features

- Reads number of equations and the equations themselves in a simple text format  
  Example:  
  `2x1 + 3x2 - x3 = 7`  
  ` -4x2 + 5x3 = -2`

- Automatically detects the highest variable index (supports x1, x2, ..., xN)

- Supports these commands/operations:

  - `num_vars` → shows the number of variables detected
  - `equation 2` → pretty-prints equation number 2
  - `column x3` → prints the entire column of coefficients for x3
  - `add 1 4` → adds equation 1 and equation 4 and shows the result
  - `solve` → solves the entire system using Gaussian elimination (only for square systems)

- Gaussian elimination with **partial pivoting** (numerically more stable than Cramer's rule)

- Handles basic coefficient parsing (including implied coefficients like `-x2` or `x5`)

## Example Input
3 

`2x1+3x2+4x3=16` 
`1x1+2x2+1x3=8`
`3x1+1x2+2x3=13`
`solve`

**Output:**

`x1 = 3`
`x2 = 2`
`x3 = 1`
