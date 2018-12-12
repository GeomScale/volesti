## SoCG 2019

- Save `liblpsolve55.so` in folder `/usr/lib`. You will find it in `/external` folder.
- In folder `/test` compile C++ code by running:  
```
cmake .  
make  
```
1. H-polytopes (Table 1):  
- To generate a unit cube in `dim` dimension and estimate the volume:  
```
./generate -cube -h -d dim
./vol -f1 cube_dim.ine -ban
```
For example:  
```
./generate -cube -h -d 20
./vol -f1 cube_20.ine -ban
```
- To generate a unit simplex in `dim` dimension and estimate the volume:  
```
./generate -simplex -h -d dim
./vol -f1 simplex_dim.ine -ban
```
- To generate a random H-polytope in `dim` with `k` facets dimension and estimate the volume:  
```
./generate -rh -d dim -m k
./vol -f1 random_h_poly_dim_k.ine -ban
```

2. V-polytopes (Table 2):  
- To generate a cross polytope in `dim` dimension and estimate the volume:  
```
./generate -cross -v -d dim
./vol -f2 cross_dim.ext -ban
```
- To generate a unit simplex in `dim` dimension:  
```
./generate -simplex -v -d dim
./vol -f2 simplex_dim.ext -ban
```
- To generate a unit cube in `dim` dimension and estimate the volume:  
```
./generate -cube -v -d dim
./vol -f2 cube_dim.ext -ban
```
- To generate a random V-polytope in `dim` dimension with k vertices and estimate the volume:  
```
./generate -rv -d dim -m k
./vol -f2 random_v_poly_dim_k.ext -ban -r
```
Note: For random V-polytopes use the flag `-r` to round the polytope.

3. Zonotopes (Table 3):  
- You can generate a random zonotope in dimension `dim` with `k` generators by running:  
```
./generate -zonotope -d dim -m k
```
- Estimate the volume using balls in MMC:  
```
./vol -f3 zonotope_dim_k.ext -ban
```
- Estimate the volume using h-polytopes in MMC:  
```
./vol -f3 zonotope_dim_k.ext -hpoly
```
- For example the following commands:  
```
./generate -zonotope -d 10 -m 15
./vol -f3 zonotope_10_15.ext -hpoly
```
Will generate a random 10-dimensional zonotope with 15 generators and estimate the volume by using h-polytopes in MMC.  
- To compute the exact volume run:  
```
./vol -f3 zonotope_10_15.ext -exact_zono
```

Note: If you wish to give a specific polytope as input use `.ine` file for an H-polytope and `.ext` file for a V-polytopes or a zonotopes. Keep the same format as in the generated files.

4. Flags

- For H-polytopes the default random walk is Coordinate Directions HnR. Use flag `-rdhr` to use Random Directions HnR:  
```
./vol -f1 cube_dim.ine -ban -rdhr
```
- For V-polytopes and zonotopes the default random walk is Random Directions HnR. To use Coordinate Directions HnR use the flag `-cdhr`. For example:  
```
./vol -f2 cube_dim.ext -ban -cdhr
```
- Use flag `-WalkL` to set the step of the HnR (the default value is 1). For example:  
```
./vol -f1 cross_dim.ext -ban -WalkL 5
```
Will set the step equals to 5.
- Use flag `-e` to set the error (the default value is `0.1`). For example:  
```
./vol -f1 zonotope_dim_k.ext -ban -e 0.2
```
- Use flag `-WinLen` to set the length of the sliding window (the default value is 4d^2+250). For example:  
```
./vol -f1 cross_dim.ext -ban -WinLen 500
```
Will set the window's length `n=500`.
- Use flags `-l` and `-u` to set the test values for testR (r) and testL (r+\delta) respectively. For example:  
```
./vol -f1 cube_dim.ine -ban -l 0.01 -u 0.015
```
Will define ratios between `0.01` and `0.015` with high probability.
- Use flag `-nuN` to set the number of points that are generated in each step of the annealing schedule, from the convex body P_i of the previous step. For example:  
```
./vol -f3 zonotope_dim_k.ext -ban -nuN 1600 10
```
Wil sample 1600 points in total and split them to 10 sub-lists. So the degrees of freedom in each t-test will be 9 = 10-1.

5. Test PCA over-aproximations of a zonotope

- To compute the ratio for the PCA approximation of a zonotope that is described in a `.ext` file, use flag `-pca` and run:  
```
./vol -f3 zonotope_dim_k.ext -hpoly -pca
```

## Volume computation and sampling

|         | Build           
| ------------- |:-------------:| 
| **master** |[![CircleCI](https://circleci.com/gh/vissarion/volume_approximation/tree/master.svg?style=svg)](https://circleci.com/gh/vissarion/volume_approximation/tree/master)
|**develop** |[![CircleCI](https://circleci.com/gh/vissarion/volume_approximation/tree/develop.svg?style=svg)](https://circleci.com/gh/vissarion/volume_approximation/tree/develop)

**VolEsti** is a C++ library for volume approximation and sampling of convex bodies (*e.g.* polytopes) with an *R* interface.

### - R Interface
------------

####  Install Rcpp package  
 
* Install package-dependencies: `Rcpp`, `RcppEigen`, `BH`, `lpSolveAPI`.  

1. Then use devtools package to install `volesti` Rcpp package. In folder /root/R-prog Run:
```r
Rcpp::compileAttributes()  
library(devtools)  
devtools::build()  
devtools::install()  
library(volesti)  
```
2. You can use Rstudio as well to open `volesti.Rproj` and then click `build source Package` and then `Install and Restart` in `Build` at the menu bar.  

#### Generate CRAN version

To generate the CRAN version of the R package follow the instructions below:  

1. From the command line navigate to folder `/cran_gen`. Then Run:  
```r
source('genCRANpkg.R')  
```

2. Open genCRANpkg.R script with `Rstudio` and run it.  

####  Run volesti from `R`
* The main function is `volume()`. It can be used to approximate the volume of a convex polytope given as a set of linear inequalities or a set of vertices (d-dimensional points) or as a Minkowski sum of segments (zonotope). There are two algorithms that can be used. The first is `SequenceOfBalls` and the second is `CoolingGaussian` (see References).  
* The function `sample_points()` can be used to sample points from a convex polytope with uniform or spherical gaussian target distribution.  
* The function `round_polytope()` can be used to round a convex polytope.  
* The function `rand_rotate()` can be used to apply a random rotation to a convex polytope.  

For more details you can read the documentation in folder `/inst/doc`.  

#### Create pdf documentation from Rd files
* Install volesti library.  
* In `R` mode (or in Rstudio) Run
```
pack = "volesti"  
path = find.package(pack)  
system(paste(shQuote(file.path(R.home("bin"), "R")),  
    "CMD", "Rd2pdf", shQuote(path)))
```
* The pdf will be created and saved in R-proj folder.  
* We give such a documentation in /R-proj/doc folder.

### - C++ Interface
------------

####  Compile C++ sources and run tests 

To compile the C++ code you have to specify the path to external library `liblpsolve55.so`, by running, in folder test:  
```
cmake -DLP_SOLVE=_PATH_TO_LIB_FILE_ .  
make  
```
For example:  `-DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so`  

You can run the tests by `cmake test` or `ctest -jK` where `K` the number of `CPU` threads. By adding the option `--verbose` to `ctest` you get more information about the tests, *e.g.* time per test, volume computed and the name of the polytope or convex body. 

#### Polytope input  

The current version of the software assumes that the polytope is given in the form of linear inequalities i.e. {x \in R^d : Ax <= b} where A is a matrix of dimension m *x* d and b a vector of dimension m or as a set of m vertices {\in R^d} or as a Minkowski sum of m segments {\in R^d}. The input is described in an `.ine`-file (H-polytopes) or in a `.ext` file (V-polytopes or zonotopes). The `.ine` file is described as follows:  
  
```  
various comments  
begin  
m d+1 numbertype  
b -A  
end  
various options  
``` 

The `.ext` file is described as follows:  
```  
various comments  
begin  
m d numbertype  
1 v_1  
.. ...  
1 v_m  
end  
various options  
``` 
In V-polytope case v_i are vertices and in zonotope case they are segments.  
  
This filestype (or similar) is used by a number of other software in polyhedral computation (e.g. `cdd`, `vinci`, `latte`). In the current version of the software, the options are treated as comments and the numbertype as C++ double type.  
If your input has equality constraints then you have to transform it in the form that only contain linear inequalities which described above by using some other software. We recommend to use latte https://www.math.ucdavis.edu/~latte for this transformation.  
  
#### Run volesti from command line  

After successful compilation you can use the software by command line. For example, the following command `./vol -h`   will display a help message about the program's available options.  
  
###### Example  
  
To estimate the volume of the 10-dimensional hypercube first prepare the file `cube10.ine` as follows:  
  
```
cube10.ine  
H-representation  
begin  
 20 11 real  
 1 1 0 0 0 0 0 0 0 0 0  
 1 0 1 0 0 0 0 0 0 0 0  
 1 0 0 1 0 0 0 0 0 0 0  
 1 0 0 0 1 0 0 0 0 0 0  
 1 0 0 0 0 1 0 0 0 0 0  
 1 0 0 0 0 0 1 0 0 0 0  
 1 0 0 0 0 0 0 1 0 0 0  
 1 0 0 0 0 0 0 0 1 0 0  
 1 0 0 0 0 0 0 0 0 1 0  
 1 0 0 0 0 0 0 0 0 0 1  
 1 -1 0 0 0 0 0 0 0 0 0  
 1 0 -1 0 0 0 0 0 0 0 0  
 1 0 0 -1 0 0 0 0 0 0 0  
 1 0 0 0 -1 0 0 0 0 0 0  
 1 0 0 0 0 -1 0 0 0 0 0  
 1 0 0 0 0 0 -1 0 0 0 0  
 1 0 0 0 0 0 0 -1 0 0 0  
 1 0 0 0 0 0 0 0 -1 0 0  
 1 0 0 0 0 0 0 0 0 -1 0  
 1 0 0 0 0 0 0 0 0 0 -1  
end  
input_incidence  
```
  
Then to use SequenceOfBalls (SOB) algorithm run the following command:  
```
./vol -f1 cube_10.ine  
```

which returns 17 numbers:  
```d m #experiments exactvolOr-1 approxVolume [.,.] #randPoints walkLength meanVol [minVol,maxVol] stdDev errorVsExact maxminDivergence time timeChebyshevBall```
  
To use CoolingGaussian (CG) algorithm run the following command:  
```
./vol -f1 cube_10.ine -CG  
```
which returns the same output as before.  

To estimate the volume of a 10-dimensional V-cross polytope described in `cross_10.ext` as follows:  
```
cross_10.ext  
V-representation  
begin  
 20 11 integer  
 1 1 0 0 0 0 0 0 0 0 0  
 1 0 1 0 0 0 0 0 0 0 0  
 1 0 0 1 0 0 0 0 0 0 0  
 1 0 0 0 1 0 0 0 0 0 0  
 1 0 0 0 0 1 0 0 0 0 0  
 1 0 0 0 0 0 1 0 0 0 0  
 1 0 0 0 0 0 0 1 0 0 0  
 1 0 0 0 0 0 0 0 1 0 0  
 1 0 0 0 0 0 0 0 0 1 0  
 1 0 0 0 0 0 0 0 0 0 1  
 1 -1 0 0 0 0 0 0 0 0 0  
 1 0 -1 0 0 0 0 0 0 0 0  
 1 0 0 -1 0 0 0 0 0 0 0  
 1 0 0 0 -1 0 0 0 0 0 0  
 1 0 0 0 0 -1 0 0 0 0 0  
 1 0 0 0 0 0 -1 0 0 0 0  
 1 0 0 0 0 0 0 -1 0 0 0  
 1 0 0 0 0 0 0 0 -1 0 0  
 1 0 0 0 0 0 0 0 0 -1 0  
 1 0 0 0 0 0 0 0 0 0 -1  
end  
hull  
incidence  
```
Run:   
```
./vol -f2 cross_10.ext  
```
which returns the same output as before.  

To estimate the volume of a 4-dimensional zonotope defined by the Minkowski sum of 8 segments described in `zonotope_4_8.ext` as follows:  
```
zonotope_4_8.ext  
Zonotpe  
begin  
 8 5 real  
 1 0.981851 -0.188734 -0.189761 0.0812645  
 1 -0.0181493 0.811266 -0.189761 0.0812645  
 1 -0.0181493 -0.188734 0.810239 0.0812645  
 1 -0.0181493 -0.188734 -0.189761 1.08126  
 1 -0.177863 0.437661 -0.0861379 -0.674634  
 1 0.737116 -0.204646 -0.540973 -0.471883  
 1 -0.684154 0.262324 0.292341 -0.265955  
 1 -0.802502 -0.740403 0.0938152 0.0874131  
end  
hull  
incidence  
```
Run:  
```
./vol -f3 zonotope_4_8.ext  
```
Flag `-v` enables the print mode.

#### Generate polytopes

You can use executable `generator` to generate polytopes (hypercubes, simplices, cross polytopes, skinny hypercubes (only in H-representation), product of two simplices (only in H-representation) and zonotoes. For example:  

1. To generate a 10-dimensional hypercube in H-representation run:  
```
./generate -cube -h -d 10
```

2. To generate a 20-dimensional simplex in V-representaion run:  
```
./generate -simplex -v -d 20
```

3. To generate a 5-dimensional zonotope defined by 10 segments run:  
```
./generate -zonotope -d 5 -m 10
```

Command `./generate -help` will display a help message about the program's available options.  

#### Sampling

You can sample from a convex polytope uniformly or from the spherical gaussian distribution. For example:  

1. To sample uniformly from the 10-dimensional hypercube, run:  
```
./vol -f1 cube_10.ine -rand -nsample 1000
```
Flag -nsample declares the number of points we wish to sample (default is 100).  

2. To sample from the gaussian distribution, run:  
```
./vol -f1 cube_10.ine -rand -nsample 1300 -gaussian -variance 1.5
```
Flag `-variance` declares the variance (default is 1.0). The center of the spherical gaussian is the Chebychev center for H-polytopes, or the origin for zonotopes. For V-polytopes is the chebychev center of the simplex that is defined by a random choice of d+1 vertices.

3. To sample from a zonotope described in zonotope.ext file run:
```
./vol -f3 zonotope.ext -rand -nsample 1500
```
For V-polytopes use flag `-f2` before the `.ext` file. In all cases use flag `-v` to print the excecutional time.  

#### Credits

Copyright (c) 2012-2018 Vissarion Fisikopoulos  
Copyright (c) 2018 Apostolos Chalkis  

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, see vol.cpp.  

The Rcpp package is based on an open-source C++ software for computing an approximation of the volume of convex bodies given as an intersection of an ellipsoid and a polytope given as an intersection of halfspaces or of a polytope given by its vertices. We have excluded CGAL and Boost dependecies. In folder include we develop the C++ code.  

Main development by Vissarion Fisikopoulos while he was affiliated with University of Athens (UoA, Greece), University of Brussels (ULB, Belgium) and Oracle Corp, and Chalkis Apostolos affiliated with University of Athens.  

#### Publications

1. I.Z. Emiris and V. Fisikopoulos, *Efficient random-walk methods for approximating polytope volume*, In Proc. ACM Symposium on Computational Geometry, Kyoto, Japan, p.318-325, 2014.  
2. I.Z. Emiris and V. Fisikopoulos, *Practical polytope volume approximation*, ACM Transactions on Mathematical Software, vol 44, issue 4, 2018.
3. L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos, *Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises*, Proc. of Symposium on Computational Geometry, Budapest, Hungary, 2018.  

