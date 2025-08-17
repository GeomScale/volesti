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
