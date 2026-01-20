## Compilation

This example demonstrates how to compute the volume of an H-polytope
using different randomized algorithms provided by volesti.

It is recommended to build the example using an out-of-source build.

```bash
mkdir build
cd build
cmake ..
make -j$(nproc)
```  
Alternatively, a local installation of lp_solve can be specified explicitly:
```bash
cmake .. -DLP_SOLVE=/path/to/liblpsolve55.so
```

## Running

After a successful build, run the executable with:
```bash
 ./hpolytopeVolume
```
The program computes the volume of a 4-dimensional H-polytope using
different randomized algorithms. Due to the randomized nature of the
methods, results may vary slightly between different runs.
