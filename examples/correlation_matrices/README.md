## Compilation
Build the example by running the following commands in this directory.

```bash
mkdir -p build && cd build
cmake .. && make
```  
HiGHS, the LP solver, is downloaded automatically by CMake.

## Usage:
```bash
./sampler && python3 ../plot.py
```