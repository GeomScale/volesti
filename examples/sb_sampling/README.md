## Explanation
In 'shake_and_bake_walk.hpp' the Running variant of Shake And Bake class of boundary sampling algorithms. It follows the steps as described in [1]. The walk can be run by using 'shakeandbake.cpp' Additionaly, to test the uniformity, the scaling ratio test ('scaling_ratio.hpp') along with random facet 2D projection for generic polytopes (cube, simplex, birkhoff) with the plots are added in 'sbtest.py' (with generated txt files). 

[1] C. G. E. Boender, R. J. Caron, J. F. McDonald, A. H. G. Rinnooy Kan,  
    H. E. Romeijn, R. L. Smith, J. Telgen i A. C. F. Vorst,  
    *Shake-And-Bake Algorithms for Generating Uniform Points on the Boundary of Bounded Polyhedra*, 1991.  
    Available at: https://doi.org/10.1016/0166-218X(91)90006-7

## Original and Limping Variant
The implemented Shake And Bake variant is Running SB because it performs significantly better than its counterparts Original and Limping. But, if one wants to test that out, here is the additional piece of code to be added in `apply` function alongside some simple enum switch / branch logic. Also, it is very important to use  `_v = GetDirection<Point>::apply(n, rng); ` (coressponds to Step 1. of Original and Limping SB from the paper ) instead of `Point v = get_direction(P,rng);` (corresponds to Step 1. for Running SB).  

```bash
NT beta;
if (mode_ == Mode::Original) {
    NT den = dot_r - dot_k;
    if (std::abs(den) < eps) continue;
    beta = std::clamp(dot_r / den, NT(0), NT(1));
} 
else {
    beta = -dot_k;
}

if (beta > NT(0) && beta <= NT(1) &&
    rng.sample_urdist() < beta)
{
    p_         = y;
    facet_idx_ = facet_new;
    A_row_k_   = A_row_r;
    Ar_.noalias() -= lambda_hit * Av_;  
}
```  
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

For sampling 10 dimensional cube with 1e-11  as epsilon (eps):
 ./billiardshakeandbake simplex 10 1e-11

 ```

 ### Output:
 ```
Parameters: walk_len=100, n_samples=2500, burn_in_iters=25 (dim=5) eps=1e-10
Generated 2500 samples in 100 steps each.
Scaling factors:
0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 

Coverage matrix (each row = one facet):
Facet 0: 0.16092 0.279693 0.386973 0.490421 0.586207 0.685824 0.762452 0.869732 0.969349 1 
Facet 1: 0.188 0.324 0.452 0.536 0.6 0.684 0.72 0.828 0.956 1   etc.


Facet        Max deviation (%)        Avg deviation (%)
     0               9.04                   6.92
     1              15.20                   7.88 etc.


In addition to diagnostics, txt file with point values across all facets will be generated as sb_[polytope]_run.
```

