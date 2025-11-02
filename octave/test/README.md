# VolEsti Octave Package Tests

This directory contains test files for the VolEsti Octave package.

## Running Tests

From within Octave, after loading the package:

```octave
pkg load volesti
runtests
```

## Individual Test Suites

- `test_gencube.m` - Tests for GenCube function
- `test_volume.m` - Tests for volume computation
- `test_sampling.m` - Tests for sampling functions
- `volesti_test.m` - Comprehensive test suite

## Test Coverage

The test suite covers:

1. **Polytope Generation**
   - Basic functionality
   - Custom scales
   - Different dimensions
   - Structure validation

2. **Volume Computation**
   - Basic volume calculation
   - Different algorithms
   - Custom error tolerances
   - Polytope structure input

3. **Sampling**
   - Basic sampling
   - Different walk methods (CDHR, RDHR, Ball, Billiard)
   - Custom walk lengths
   - Burn-in steps

4. **Error Handling**
   - Invalid inputs
   - Dimension mismatches
   - Missing fields

## Automated Testing

Tests can be run automatically as part of CI/CD pipelines. Example:

```bash
octave --eval "pkg load volesti; runtests"
```

