## Compilation
Build the example by running the following commands in this directory.

```bash
mkdir -p build && cd build
cmake . -DLP_SOLVE=_PATH_TO_LPSOLVE_LIB && make
```  
You have to specify the path to liblpsolve55.so/dll/dylib.  
For example: -DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so

## Usage:
```bash
 build/sampler
 python3 plot.py
```