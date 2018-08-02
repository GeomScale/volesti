## Volume computation and sampling

|         | Build           
| ------------- |:-------------:| 
| **master** |[![CircleCI](https://circleci.com/gh/vissarion/volume_approximation/tree/master.svg?style=svg)](https://circleci.com/gh/vissarion/volume_approximation/tree/master)
|**develop** |[![CircleCI](https://circleci.com/gh/vissarion/volume_approximation/tree/develop.svg?style=svg)](https://circleci.com/gh/vissarion/volume_approximation/tree/develop)

**VolEsti** is a C++ library for volume approximation and sampling of convex bodies (*e.g.* polytopes) with an *R* interface.

Documentation
----------------

####  Install Rcpp package  
 
* Install package-dependencies: `Rcpp`, `RcppEigen`, `BH`, `lpSolveAPI`.  
* Then use devtools package to install `volesti` Rcpp package. Run:
```r
Rcpp::compileAttributes()  
build()  
install()  
library(volesti)  
```
 You can use Rstudio as well to open `volesti.Rproj` and then click `build source Package` and then `Install and Restart` in `Build` at the menu bar.  

####  Run volesti from `R`
* The main function is `VolEsti()`. The input has to be a `list("matrix"=A, "vector"=b, "cheb"=xc, "rounding"=bool, "verbose"=BOOL)`, for a polytope Ax<=b and a d+1 vector `xc` which last coordinate is the radius of the chebychev ball and the first d coordinates the center. Rounding option is to apply a linear transformation to the convex body to get a well rounded one.  
* You can give as input a `list("path"='path/to/ine/file')` insteed of a `"matrix"` or a `"vector"`.  
* The `"cheb"` input is optional. When it is not given lpsolve library is used from C++ code.  
* An `R` function, using `lpSolveAPI`, to compute Chebychev center is provided.
* `"verbose"` is by default `false`.  
* You can run all the tests by running function `testRvolEsti()`.  

#### Create pdf documentation from Rd files
* Navigate to /R-proj folder from the commend line.  
* Run
```
R CMD Rd2pdf --pdf --title='VolEsti Documentation' -o /path/to/save/volesti.pdf man/*.Rd
```
* The pdf will be created and saved in the declared path.  
* We give such a documentation in /R-proj/doc folder.

####  Compile C++ sources and run tests 

To compile the C++ code you have to specify the path to external library `liblpsolve55.so`, by running, in folder test:  
```
cmake -DLP_SOLVE=_PATH_TO_LIB_FILE_ .  
make  
```
For example:  `-DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so`  

You can run the tests by `cmake test` or `ctest -jK` where `K` the number of `CPU` threads. By adding the option `--verbose` to `ctest` you get more information about the tests, *e.g.* time per test, volume computed and the name of the polytope or convex body. 

#### Polytope input  

The current version of the software assumes that the polytope is given in the form of linear inequalities i.e. {x \in R^d : Ax <= b} where A is a matrix of dimension m *x* d and b a vector of dimension m. The input is described in an `.ine`-file as follows:  
  
```  
various comments  
begin  
m d+1 numbertype  
b -A  
end  
various options  
``` 
  
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
  
Then to use volesti algorithm run the following command:  
```
./vol -f1 ./test_data/cube10.ine  
```

which returns 17 numbers:  
```d m #experiments exactvolOr-1 approxVolume [.,.] #randPoints walkLength meanVol [minVol,maxVol] stdDev errorVsExact maxminDivergence time timeChebyshevBall```
  
To use CV algorithm run the following command:  
```
./vol -f1 ./test_data/cube10.ine -g_an  
```
which returns:  
```volume computed =...  
Total time = .. sec```
 
Credits
-------

Copyright (c) 2012-2018 Vissarion Fisikopoulos  
Copyright (c) 2018 Apostolos Chalkis

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, see vol.cpp.  

The Rcpp package is based on an open-source C++ software for computing an approximation of the volume of convex bodies given as an intersection of an ellipsoid and a polytope given as an intersection of halfspaces or of a polytope given by its vertices. We have excluded CGAL and Boost dependecies. In folder include we develop the C++ code.

Main development by Vissarion Fisikopoulos while he was affiliated with University of Athens (UoA, Greece), University of Brussels (ULB, Belgium) and Oracle Corp, and Chalkis Apostolos affiliated with University of Athens.

#### Publications

1. I.Z. Emiris and V. Fisikopoulos, *Efficient random-walk methods for approximating polytope volume*, In Proc. ACM Symposium on Computational Geometry, Kyoto, Japan, p.318-325, 2014.  
2. I.Z. Emiris and V. Fisikopoulos, *Practical polytope volume approximation*, ACM Transactions on Mathematical Software, vol 44, issue 4, 2018.
3. L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos, *Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises*, Proc. of Symposium on Computational Geometry, Budapest, Hungary, 2018.


