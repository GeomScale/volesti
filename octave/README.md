# VolEsti Octave Package

This directory contains the Octave interface for the VolEsti library.

## Overview

VolEsti is a C++ library for volume approximation and sampling of convex bodies (e.g. polytopes). This package provides an Octave interface that allows Octave users to:

- Compute volumes of convex polytopes
- Sample uniform points from polytopes
- Generate common polytope types (e.g., hypercubes)

## Installation

### Prerequisites

- Octave (version 5.1.0 or later)
- C++ compiler with C++17 support
- CMake (>= 3.11)
- Eigen3
- Boost libraries
- lp_solve

### Building MEX Files

1. Navigate to the source directory:
   ```bash
   cd octave/src
   mkdir build
   cd build
   ```

2. Configure and build:
   ```bash
   cmake ..
   make
   ```

   On Windows with MSVC:
   ```cmd
   cmake ..
   cmake --build . --config Release
   ```

3. The compiled `.oct` files should be in `octave/inst/`.

### Installing the Package

From within Octave, navigate to the `octave` directory and run:
```octave
pkg install .
```

To load the package:
```octave
pkg load volesti
```

## Usage Examples

### Computing Volume

```octave
% Generate a 3-dimensional cube
P = GenCube(3);
% Compute its volume (should be close to 8)
vol = volume(P);
printf("Volume: %f\n", vol);
```

### Sampling Points

```octave
% Generate a cube
P = GenCube(3);
% Sample 100 uniform points
points = sample_points(P, 100);
% Plot the first two dimensions
plot(points(1,:), points(2,:), 'o');
```

### Custom Polytope

```octave
% Define a polytope: -1 <= x_i <= 1 for i=1,2,3
A = [eye(3); -eye(3)];
b = ones(6, 1);
% Compute volume
vol = volume(A, b);
% Sample points
points = sample_points(A, b, 50);
```

## Package Structure

- `DESCRIPTION` - Package metadata
- `COPYING` - License information
- `inst/` - Installed files (Octave functions)
- `src/` - Source code (MEX files)
- `doc/` - Documentation

## Functions

- `volume()` - Compute volume of a polytope
- `sample_points()` - Sample uniform points from a polytope
- `GenCube()` - Generate a hypercube polytope

## Testing

The package includes a comprehensive test suite. To run all tests:

```octave
pkg load volesti
runtests
```

You can also run individual test suites:
```octave
test_gencube      % Test polytope generation
test_volume       % Test volume computation
test_sampling     % Test sampling functions
volesti_test      % Comprehensive test suite
```

## Documentation

See the Octave help system:
```octave
help volume
help sample_points
help GenCube
```

## License

This package is distributed under the GNU Lesser General Public License, version 3 or later.
See the LICENSE file in the parent directory.

## Contributing

Contributions are welcome! Please see the main CONTRIBUTING.md file for guidelines.

## References

- [VolEsti GitHub Repository](https://github.com/GeomScale/volesti)
- [VolEsti Documentation](https://volesti.readthedocs.io)
- [GeomScale Project](https://geomscale.github.io)

