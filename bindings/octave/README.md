# GNU Octave Interface for Volesti - Proof of Concept

This directory contains a **Proof of Concept** implementation of a native GNU Octave interface for the [Volesti](https://github.com/GeomScale/volesti) high-performance computational geometry library.

## Overview

Octave bindings for Volesti's computational geometry library, following the same API design as [Rvolesti](https://github.com/GeomScale/Rvolesti) for consistency across R and Octave.

### Features

- **Volume Computation**: High-dimensional polytope volume using Sequence of Balls (SOB) algorithm
- **H-Polytope Support**: Polytopes defined by Ax ≤ b (half-space representation)
- **V-Polytope Support**: Polytopes defined by convex hull of vertices  
- **Memory Efficient**: The interface uses `Eigen::Map` to minimize memory copying.
- **Rvolesti-Compatible API**: Same constructor and function names as R bindings

### Directory Structure

```
bindings/octave/
├── src/          # C++ MEX source files
├── inst/         # User-facing API (M-files)
├── test/         # Test scripts
└── examples/     # Demo scripts
```


## Prerequisites

### Required Software

- **GNU Octave** (≥ 4.0) with development headers
- **C++ Compiler** with C++17 support (g++ or clang++)
- **Eigen** library (included in Volesti's `external/` directory)
- **Boost** library (included in Volesti's `external/` directory)

### Installation

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install -y octave liboctave-dev
```

#### Fedora/RHEL
```bash
sudo dnf install -y octave octave-devel
```

#### Arch Linux
```bash
sudo pacman -S octave
```

#### Verify Installation
```bash
octave --version
mkoctfile --version
```

## Build Instructions

### Prerequisite: Populate External Dependencies

Run CMake from the repository root to download dependencies:

```bash
cd /path/to/volesti
cmake -S . -B build
```

This populates `external/_deps/` with Eigen and Boost.

### Building the Bindings

```bash
cd bindings/octave
make
```

This compiles `src/volume.cpp` → `inst/volume.oct`.

### Running Tests

```bash
# Run all tests (like ctest)
make test

# Or run individual tests
octave --eval "addpath('inst'); cd('test'); test_Hvol"
octave --eval "addpath('inst'); cd('test'); test_Vvol"
```

## Usage

### API Functions

#### `Hpolytope(A, b)` - Create H-Polytope

Creates an H-polytope struct representing Ax ≤ b.

**Example:**
```matlab
addpath('inst');

% Create a 2D unit square [-1,1]²
A = [1 0; -1 0; 0 1; 0 -1];
b = ones(4, 1);
P = Hpolytope(A, b);
```

#### `Vpolytope(V)` - Create V-Polytope

Creates a V-polytope struct from vertices (convex hull).

**Example:**
```matlab
% Create a 3D cube from 8 vertices
V = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1;
     -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
P = Vpolytope(V);
```

#### `volesti_volume(P [, epsilon, walk_length, verbose])` - Compute Volume

Computes polytope volume using Sequence of Balls algorithm.

**Parameters:**
- `P` : Polytope struct (from `Hpolytope()` or `Vpolytope()`)
- `epsilon` : (optional) Error tolerance (default: 1.0, smaller = more accurate)
- `walk_length` : (optional) Random walk length (default: 1, larger = better mixing)
- `verbose` : (optional) Show progress (default: true)

**Examples:**
```matlab
addpath('inst');

% H-polytope: 3D cube
A = [eye(3); -eye(3)];
b = ones(6, 1);
P = Hpolytope(A, b);
vol = volesti_volume(P);  % ~8.0

% V-polytope: 3D cube from vertices  
V = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1;
     -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
P = Vpolytope(V);
vol = volesti_volume(P);  % ~8.0

% High accuracy
vol = volesti_volume(P, 0.1, 10);
```

**Note:** Function is named `volesti_volume()` instead of `volume()` to avoid conflict with Octave's built-in `volume()` function.

### Polytope Generators

#### `gen_cube(d, type)` - Generate Hypercube

Creates a d-dimensional unit hypercube [-1,1]^d.

```matlab
P = gen_cube(10, 'H');  % 10D cube in H-representation
P = gen_cube(5, 'V');   % 5D cube in V-representation
vol = volesti_volume(P);  % Returns ~2^d
```

#### `gen_cross(d, type)` - Generate Cross Polytope

Creates a d-dimensional cross polytope (dual of hypercube).

```matlab
P = gen_cross(5, 'H');   % 5D cross in H-representation
P = gen_cross(15, 'V');  % 15D cross in V-representation
vol = volesti_volume(P);  % Returns ~2^d / d!
```

#### `gen_prod_simplex(d)` - Generate Product of Simplices

Creates a 2d-dimensional polytope as product of two d-dimensional simplices (H-representation only).

```matlab
P = gen_prod_simplex(5);  % Product of two 5D simplices (10D polytope)
vol = volesti_volume(P);  % Returns ~(1/d!)^2
```

## Technical Architecture

### Single-Copy Memory Architecture

Traditional language bindings copy data between the host language and the library, which is expensive for large matrices. The interface uses `Eigen::Map` to minimize memory copying:

**Traditional approach (2 copies):**
```cpp
// Octave data -> C++ temporary (copy 1)
std::vector<double> A_cpp(octave_A.data(), octave_A.data() + octave_A.numel());
// C++ temporary -> Eigen (copy 2)
Eigen::MatrixXd A_eigen = Eigen::Map<Eigen::MatrixXd>(A_cpp.data(), m, n);
```

**Our approach (1 copy total):**
```cpp
Eigen::Map<MT> A_eigen(octave_A.fortran_vec(), m, n);  // Zero-copy view into Octave memory
Hpolytope P(n, A_eigen, b_eigen);  // Single copy into polytope's internal storage
```

**Benefits:**
- **O(1) time complexity** for data transfer (vs O(n) for copying)
- **~50% reduction** in memory allocation overhead
- **Cache-friendly** direct memory access

### Build System

The Makefile uses `mkoctfile`, Octave's wrapper around g++/clang++ that handles:
- Octave header paths
- Shared library linking
- Platform-specific compilation flags

## Known Limitations

- **Unbounded polytopes**: Not detected, will fail
- **Empty polytopes**: May crash during interior point computation  
- **High dimensions**: Performance degrades significantly above ~50D
- **Numerical precision**: Very small (<1e-8) or large (>1e8) polytopes may have issues
- **Sparse matrices**: Not optimized (uses dense storage)
- **Single algorithm**: Only exposes `volume_sequence_of_balls`

These limitations will be addressed in the full implementation.

## License

This code follows Volesti's license (GNU LGPL 3.0).

## Contact

For questions or issues related to this PoC, please open an issue on the [Volesti GitHub repository](https://github.com/GeomScale/volesti).

