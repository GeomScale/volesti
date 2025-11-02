# Pull Request: GNU Octave Interface for VolEsti

## Summary

This PR adds a complete GNU Octave package interface for VolEsti, allowing Octave users to use the library with native Octave code.

## Changes

### Package Structure
- Complete Octave package structure following Octave package conventions
- `DESCRIPTION`, `COPYING`, `INDEX` files for package metadata
- Installation instructions (`INSTALL`, `README.md`)

### Implementation
- **MEX Interfaces** (`src/`):
  - `volesti_volume_mex.cpp` - Volume computation MEX interface
  - `volesti_sample_mex.cpp` - Sampling MEX interface
  - `CMakeLists.txt` - Build configuration for MEX files
  - Build scripts for Linux/Windows

- **Octave Functions** (`inst/`):
  - `volume()` - Compute volume of polytopes
  - `sample_points()` - Sample uniform points from polytopes
  - `GenCube()` - Generate hypercube polytopes
  - `volesti()` - Package information function
  - `volesti_example()` - Example script

- **Test Suite** (`inst/`):
  - `volesti_test()` - Comprehensive test suite (11 test cases)
  - `test_gencube()` - Polytope generation tests
  - `test_volume()` - Volume computation tests
  - `test_sampling()` - Sampling tests
  - `runtests()` - Test runner

### Documentation
- Texinfo documentation (`doc/volesti.texi`)
- Updated main `README.md` with Octave interface information
- Updated `docs/getting_started/install.md` with installation instructions
- Package README with usage examples

## Features

### Volume Computation
- Multiple algorithms: sequence_of_balls, cooling_gaussians
- Customizable error tolerance
- Support for H-polytopes (Ax <= b)

### Sampling
- Multiple random walk methods:
  - CDHR (Coordinate Directions Hit-and-Run)
  - RDHR (Random Directions Hit-and-Run)
  - Ball Walk
  - Billiard Walk
- Customizable walk length and burn-in steps

### Polytope Generation
- Hypercube generator with customizable dimensions and scale

## Testing

Comprehensive test suite included:
- Functionality tests for all main functions
- Error handling and input validation
- Edge cases (different dimensions, parameters)
- All tests can be run with `runtests()`

## Installation

See `octave/README.md` and `octave/INSTALL` for detailed installation instructions.

## Status

✅ Package structure complete
✅ MEX interfaces implemented
✅ Octave wrapper functions implemented
✅ Test suite created
✅ Documentation added
✅ Build system configured
⚠️ Needs compilation and testing on target platforms

## Next Steps

1. Build and test on Linux, macOS, and Windows
2. Verify compatibility with different Octave versions
3. Submit to Octave Packages index after testing
4. Add CI/CD for automated testing

