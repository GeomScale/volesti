# SuiteSparse Integration

This directory contains CMake files for downloading and building external dependencies, including SuiteSparse.

## SuiteSparse

The `SuiteSparse.cmake` file handles the automatic download, build, and configuration of SuiteSparse from the official GitHub repository.

### Components Built

The following SuiteSparse components are built:
- **SuiteSparse_config**: Core configuration
- **AMD**: Approximate Minimum Degree ordering
- **COLAMD**: Column Approximate Minimum Degree ordering  
- **CAMD**: Constrained Approximate Minimum Degree ordering
- **CCOLAMD**: Constrained Column Approximate Minimum Degree ordering
- **CHOLMOD**: Sparse Cholesky factorization
- **SPQR**: Sparse QR factorization

### CMake Integration

The SuiteSparse library is fetched and built automatically when you configure the project:

```bash
cd test
mkdir build
cd build
cmake ..
make
```

### Dependencies

SuiteSparse requires:
- BLAS (Basic Linear Algebra Subprograms)
- LAPACK (Linear Algebra Package)

These are typically available on most systems:
- **macOS**: Pre-installed as part of Accelerate framework
- **Linux**: Install via `libblas-dev` and `liblapack-dev` packages

The CMake build system will automatically find these dependencies.

### Version

The integration uses SuiteSparse release: `v7.12.1`

This ensures a stable, tested version of the library.

### Build Output

After building, SuiteSparse libraries are installed to:
```
external/_deps/suitesparse-build/install/
├── include/
│   └── suitesparse/
└── lib/
    ├── libspqr.a
    ├── libcholmod.a
    ├── libamd.a
    ├── libcolamd.a
    ├── libcamd.a
    ├── libccolamd.a
    └── libsuitesparseconfig.a
```
