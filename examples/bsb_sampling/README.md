## Billiard Shake and Bake for H Polytopes
This is an example that illustrates the simplification and preparation process.

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
 ./billiardshakeandbake <cube|simplex|birkhoff> <dimension> [nr] [epsilon]
```

## Example:
```
For sampling 10 dimensional cube with default upper bound for reflections (nr) and epsilon (eps):
 ./billiardshakeandbake cube 10

For sampling 10 dimensional cube with 2 as upper bound for reflections (nr) and default epsilon (eps):
 ./billiardshakeandbake simplex 10 2

For sampling 10 dimensional cube with default upper bound for reflections (nr) and 1e-11 as epsilon (eps):
 ./billiardshakeandbake simplex 10 -1 1e-11

 ```

 ### Output:
 ```
Parameters: walk_len=500, n_samples=12500, burn_in_iters=125 (dim=25, nr=5) eps=1e-10
Generated 12500 samples in 500 steps each.
Scaling factors:
0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 

Coverage matrix (each row = one facet):
Facet 0: 0.127413 0.216216 0.34749 0.447876 0.548263 0.667954 0.764479 0.837838 0.915058 1 
Facet 1: 0.134387 0.221344 0.3083 0.407115 0.490119 0.581028 0.695652 0.794466 0.909091 1  etc.


Facet        Max deviation (%)        Avg deviation (%)
     0               6.80                   3.73
     1               3.44                   1.19 etc.


In addition to diagnostics, txt file with point values across all facets will be generated as billiard_sb_[polytope]_[nr].
```

