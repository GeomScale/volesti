# Full Dimensional Polytope Examples

This directory contains examples demonstrating the `compute_full_dimensional_polytope` function from `include/preprocess/full_dimensional_polytope.hpp`.

## Overview

The `compute_full_dimensional_polytope` function transforms a potentially lower-dimensional polytope defined by:
- **Equality constraints**: `Aeq * x = beq`
- **Inequality constraints**: `A * x <= b`

Into a full-dimensional polytope: `A_full * y <= b_full`

The transformation provides:
- **shift**: A point satisfying the equality constraints
- **N**: Nullspace basis matrix with orthonormal columns
- **Mapping**: `x = shift + N * y` (from reduced space to original space)

## Building the Example

From this directory:

```bash
mkdir build
cd build
cmake ..
make
```

## Running the Example

```bash
./full_dimensional_example
```

The program will run through 6 examples, pausing between each for you to review the output.

## Output Format

For each example, the program displays:

1. **Description**: What the example demonstrates
2. **Input**: The equality and inequality constraints
3. **Results**:
   - Original and reduced dimensions
   - Shift vector
   - Nullspace basis N
   - Output polytope dimensions
   - Orthonormality verification

## Key Features Demonstrated

- ✅ Automatic dimension reduction
- ✅ Orthonormal nullspace basis
- ✅ Handling of sparse constraints
- ✅ Multiple equality constraints
- ✅ Extreme reductions (nD → 1D)
- ✅ Real-world polytopes (Birkhoff)

## Mathematical Background

Given a polytope `P = {x : Aeq*x = beq, A*x <= b}` which is lower-dimensional (rank(Aeq) > 0), the function computes:

1. **Shift vector**: Solves `Aeq * shift = beq` using SPQR
2. **Nullspace basis**: Computes `N` such that `Aeq * N = 0` via QR factorization
3. **Transformed polytope**: `P' = {y : A*N*y <= b - A*shift}`

Any point `y` in `P'` maps to a feasible point in `P` via: `x = shift + N*y`

## Dependencies

- Eigen3 (matrix operations)
- SuiteSparse (SPQR for QR factorization, CHOLMOD for linear systems)
- Boost (for polytope generators)
- lp_solve (for polytope operations)

## See Also

- **Header file**: `include/preprocess/full_dimensional_polytope.hpp`
- **Tests**: `test/full_dimensional_polytope_test.cpp`
