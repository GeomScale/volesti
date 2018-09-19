# Volume computation and sampling
----------------------------------

###  About  
[**volesti**](https://github.com/vissarion/volume_approximation) is a C++ library for volume approximation and sampling of convex bodies (*e.g.* polytopes). This package provides a [CGAL](https://www.cgal.org/) free version and a *R* interface.  

Package `volesti` computes approximations of volume of polytopes given as a set of points or linear inequalities or as a Minkowski sum of segments (zonotopes). There are two algorithms for volume approximation as well as algorithms for sampling, rounding and rotating polytopes.  

###  Install CRAN package
* Install `volesti` by runing:  
```
install.packages("volesti")
```
* The package-dependencies are: `Rcpp`, `RcppEigen`, `BH`.  

###  Usage
* The main function is `volume()`. It can be used to approximate the volume of a convex polytope given as a set of linear inequalities or a set of vertices (d-dimensional points) or as a Minkowski sum of segments (zonotope). There are two algorithms that can be used. The first is `SequenceOfBalls` and the second is `CoolingGaussian` (see References).  
* The function `sample_points()` can be used to sample points from a convex polytope with uniform or spherical gaussian target distribution.  
* The function `round_polytope()` can be used to round a convex polytope.  
* The function `rand_rotate()` can be used to apply a random rotation to a convex polytope.  

For more details you can read the documentation in folder `/inst/doc`.  

### Credits
-------

Copyright (c) 2012-2018 Vissarion Fisikopoulos  
Copyright (c) 2018 Apostolos Chalkis  

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  

In folder `src/include` we develop the C++ code.

Main development by Vissarion Fisikopoulos while he was affiliated with University of Athens (UoA, Greece), University of Brussels (ULB, Belgium) and Oracle Corp, and Chalkis Apostolos affiliated with University of Athens.

### Publications

1. I.Z. Emiris and V. Fisikopoulos, *Efficient random-walk methods for approximating polytope volume*, In Proc. ACM Symposium on Computational Geometry, Kyoto, Japan, p.318-325, 2014.  
2. I.Z. Emiris and V. Fisikopoulos, *Practical polytope volume approximation*, ACM Transactions on Mathematical Software, vol 44, issue 4, 2018.  
3. L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos, *Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises*, Proc. of Symposium on Computational Geometry, Budapest, Hungary, 2018.  

