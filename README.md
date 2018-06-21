## *VolEsti:* Approximate computation of volume of convex bodies.

Copyright (c) 2012-2017 Vissarion Fisikopoulos

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, see vol.cpp.

Main development by Vissarion Fisikopoulos, University of Athens, Greece.
Main algorithm based on: I.Z. Emiris and V. Fisikopoulos, "Efficient random-walk methods for approximating polytope volume", In Proc. ACM Symposium on Computational Geometry 2014, pages 318-325.

### *We develop the R version of VolEsti*  

* We develop a C++ library and a R interface.
* We have excluded CGAL and Boost dependecies.
* In folder include we develop the C++ code.
* To run C++ code you have to give the path to external library liblpsolve55.so in the CMakeLists.txt file.
* Then  run:  
cmake .  
make
* To run the R interface you have to run RcppExports.R script, then the command Rcpp::sourceCpp('path/to/R-proj/src/vol_R.cpp'). Then VolEsti() R function can be used by giving a list("matrix"=A, "vector"=b, "cheb"=xc), for a polytope Ax<=b and a d+1 vector xc which last coordinate is the radius of the chebychev ball and the first d coordinates the center.
* The "cheb" input is optional. When it is not given lpsolve is used.
* The "vector" input is optional when A is at ".ine" format and the first row is (m,d,0, ... ,0) where m is the number of facets and d the dimension.
* In folder R-proj/src we develop the R interface. You can call the C++ library through R using command: Rcpp::sourceCpp('path/to/R-proj/src/vol_R.cpp'). Then function vol_R() is ready to call.
* You have to install Rcpp, RcppEigen and lpSolveAPI R libraries.

