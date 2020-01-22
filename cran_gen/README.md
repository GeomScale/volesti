# Volume computation and sampling

## About  
The `volesti` package provides [R](https://www.r-project.org/) with functions for volume estimation and sampling. In particular, it provides an R interface for the C++ library [**volesti**](https://github.com/GeomScale/volume_approximation). 

`volesti` computes approximations of volume of polytopes given as a set of points or linear inequalities or as a Minkowski sum of segments (zonotopes). There are two algorithms for volume approximation as well as algorithms for sampling, rounding and rotating polytopes.  

##  Download and install 

* The latest stable version is available from CRAN.
* The latest development version is available on Github `www.github.com/GeomScale/volume_approximation`

* Install `volesti` by running:  
```
install.packages("volesti")
```
* The package-dependencies are: `Rcpp`, `RcppEigen`, `BH`. 

##  Usage

* The main function is `volume()`. It can be used to approximate the volume of a convex polytope given as a set of linear inequalities or a set of vertices (d-dimensional points) or as a Minkowski sum of segments (zonotope). There are two algorithms that can be used. The first is `SequenceOfBalls` and the second is `CoolingGaussian`.  
* `sample_points()` can be used to sample points from a convex polytope approximating uniform or multidimensional spherical gaussian target distribution using: (a) coordinate directions hit-and-run (default), (b) random directions hit-and-run or (c) ball walk. This function can be used as well in order to sample exact uniform points from simplices and the boundary or the interior of hyperspheres.  
* `round_polytope()` can be used to round a convex polytope.  
* `rand_rotate()` can be used to apply a random rotation to a convex polytope.  
* `gen_cross()` can be used to generate a d-dimensional cross polytope.  
* `gen_cube()` can be used to generate a d-dimensional unit cube.  
* `gen_prod_simplex()` can be used to generate a 2d-dimensional polytope that is defined as the product of two d-dimensional unit simplices.  
* `gen_simplex()` can be used to generate the d-dimensional unit simplex.  
* `gen_skinny_cube()` can be used to generate a d-dimensional skinny hypercube.  
* `gen_rand_hpoly()` can be used to generate a d-dimensional polytope in H-representation with m facets.  
* `gen_rand_vpoly()` can be used to generate a d-dimensional polytope in V-representation with m vertices.  
* `gen_rand_zonotope()` can be used to generate a random d-dimensional zonotope defined by the Minkowski sum of m d-dimensional segments.  
* `file_to_polytope()` takes the path for an ine or an ext file and returns the corresponding polytope.  
* `frustum_of_simplex()` can be used to evaluate the portion of of the d-dimensional unit simplex contained in a given half-space using M. Ali's version of G. Varsi's iterative formula.  
* `inner_ball()` can be used to compute an inscribed ball of a convex polytope (works for all the representations).  
* `copula1()` can be used to compute the copula when two families of parallel hyperplanes are given.  
* `copula2()` can be used to compute the copula when one familiy of parallel hyperplanes  and a family of concentric ellipsoids are given.  

For more details, features, examples and references you can read the documentation.  

## Credits

Copyright (c) 2012-2020 Vissarion Fisikopoulos  
Copyright (c) 2018-2020 Apostolos Chalkis  

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.   

Main development by Vissarion Fisikopoulos while he was affiliated with University of Athens (UoA, Greece) and University of Brussels (ULB, Belgium), and Chalkis Apostolos affiliated with University of Athens. Part of the development was done while  A.Chalkis (as student) and V.Fisikopoulos (as mentor) were participating in Google Summer of Code 2018 program.

### Publications

1. I.Z. Emiris and V. Fisikopoulos, *Efficient random-walk methods for approximating polytope volume*, In Proc. ACM Symposium on Computational Geometry, Kyoto, Japan, p.318-325, 2014.  
2. I.Z. Emiris and V. Fisikopoulos, *Practical polytope volume approximation*, ACM Transactions on Mathematical Software, vol 44, issue 4, 2018.  
3. L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos, *Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises*, Proc. of Symposium on Computational Geometry, Budapest, Hungary, 2018.  

