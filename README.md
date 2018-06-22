## *VolEsti:* Approximate computation of volume of convex bodies.

Copyright (c) 2012-2017 Vissarion Fisikopoulos

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, see vol.cpp.

Main development by Vissarion Fisikopoulos, University of Athens, Greece.
Main algorithm based on: I.Z. Emiris and V. Fisikopoulos, "Efficient random-walk methods for approximating polytope volume", In Proc. ACM Symposium on Computational Geometry 2014, pages 318-325.

### *We develop the R version of VolEsti*  

* We develop a C++ library and a R interface.
* We have excluded CGAL and Boost dependecies.
* In folder include we develop the C++ code.
* To run C++ code you have to specify the path to external library liblpsolve55.so, by running:  
cmake -DLP_SOLVE=_PATH_TO_LIB_FILE_ .  
make
* For example:  
cmake /usr/lib/lpsolve/liblpsolve55.so .  
make  

* To use R interface you have to install Rcpp, RcppEigen and lpSolveAPI R libraries.
* To run the R interface (turn on library('lpSolveAPI')) you have two choices:  
1. With Rstudio: Open the project in folder R-proj and click "Build source package" from Build menu at the toolbar, and then click "Install and restart" from Build menu as well.  
2. From command line: Install R library "devtools". Navigate to /path/to/R-proj and run:  
>$R  
>Rcpp::compileAttributes()  
>build()  
>install()  
>library(volesti)  
* The main function is VolEsti. The input has to be a list("matrix"=A, "vector"=b, "cheb"=xc), for a polytope Ax<=b and a d+1 vector xc which last coordinate is the radius of the chebychev ball and the first d coordinates the center.
* The "cheb" input is optional. When it is not given lpsolve library is used.  
* "vector"=b is optional as well: You can also give as an input only a matrix A, when it is in the same format with matrix "A3.Rdata" in src folder (".ine" style).

