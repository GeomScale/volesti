## *VolEsti:* Approximate computation of volume of convex bodies.

Copyright (c) 2018 Vissarion Fisikopoulos, Chalkis Apostolos

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, see vol.cpp.  

Main development by Vissarion Fisikopoulos, Oracle Greece, and Chalkis Apostolos, University of Athens, Greece.
Main algorithm based on: I.Z. Emiris and V. Fisikopoulos, "Efficient random-walk methods for approximating polytope volume", In Proc. ACM Symposium on Computational Geometry 2014, pages 318-325.  

This is a Rcpp package based on an open-source C++ software for computing an approximation of the volume of convex bodies given as an intersection of an ellipsoid and a polytope given as an intersection of halfspaces or of a polytope given by its vertices. We have excluded CGAL and Boost dependecies. In folder include we develop the C++ code.

### *Install Rcpp package*  

* Install package-dependencies: Rcpp, RcppEigen, BH, lpSolveAPI.  
* Then use devtools package to install volesti Rcpp package. Run:
>>Rcpp::compileAttributes()  
>build()  
>install()  
>library(volesti)  
You can use Rstudio as well to open volesti.Rproj and then click "build source Package" and then "Install and Restart" in "Build" at the menu bar.  

#### *Run volesti*

* The main function is VolEsti(). The input has to be a list("matrix"=A, "vector"=b, "cheb"=xc, "verbose"=BOOL), for a polytope Ax<=b and a d+1 vector xc which last coordinate is the radius of the chebychev ball and the first d coordinates the center.  
* You can give as input a list("path"='path/to/ine/file') insteed of a "matrix" or a "vector".  
* The "cheb" input is optional. When it is not given lpsolve library is used from C++ code.  
* A R function, using lpSolveAPI, to compute chebychev center is provided.
* Verbose variable is by default false.  
* You can run all the experiments by running function testRvolEsti().  

### *Compile C++ Sources*  

* To run C++ code you have to specify the path to external library liblpsolve55.so, by running, in folder test:  
>cmake -DLP_SOLVE=_PATH_TO_LIB_FILE_ .  
>make  
* For example:  
cmake -DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so .  
make  
* You can compile the c++ code directly by running ./run_tests.sh. If you give input -exp then after compile all tests in test_data folder will be runned and declared as PASSED or FAILLED.  
  
#### *Polytope input*  

* The current version of the software assumes that the polytope is given in the form of linear inequalities i.e. {x \in R^d : Ax <= b} where A is a matrix of dimension mxd and b a vector of dimension m. The input is described in an .ine-file as follows:  
  
===================   
various comments  
begin  
m d+1 numbertype  
b -A  
end  
various options  
===================  
  
* This filestype (or similar) is used by a number of other software in polyhedral computation (e.g. cdd, vinci, latte). In the current version of the software, the options are treated as comments and the numbertype as C++ double type.  
* If your input has equality constraints then you have to transform it in the form that only contain linear inequalities which described above by using some other software. We recommend to use latte (https://www.math.ucdavis.edu/~latte) for this transformation.  
  
#### *Use volume approximation*  

After successful compilation you can use the software by command line. 
 
./vol -h  
 
will display a help message about the program's available options.  
  
**Example**  
  
* To estimate the volume of the 10-dimensional hypercube first prepare the file cube10.ine as follows:  
  
======================  
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
=======================  
  
* Then run the following command:  
>./vol -f1 polytope_examples/cube10.ine  
which returns 17 numbers:  
>d m #experiments exactvolOr-1 approx [.,.] #randPoints walkLength meanVol [minVol,maxVol] stdDev errorVsExact maxminDivergence time timeChebyshevBall
 



