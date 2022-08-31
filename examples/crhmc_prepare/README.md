## Constrained Riemannian Hamiltonian Monte Carlo for Polytopes
This is an example that illustrates the simplification and preparation process.

References:
Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with Riemannian Hamiltonian
 Monte Carlo in a Constrained Space"
## Compilation
Build the example by running the following commands in this directory.

```bash
cmake . -DLP_SOLVE=_PATH_TO_LIB_FILE
make
```  
You have to specify the path to liblpsolve55.so/dll/dylib.  
For example: -DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so

## Running:
```bash
 ./crhmc_prepare [filepath]
```
