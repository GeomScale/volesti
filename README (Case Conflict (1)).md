## *VolEsti:* Approximate computation of volume of convex bodies.

Copyright (c) 2012-2017 Vissarion Fisikopoulos

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, see vol.cpp.

Main development by Vissarion Fisikopoulos, University of Athens, Greece.
Main algorithm based on: I.Z. Emiris and V. Fisikopoulos, "Efficient random-walk methods for approximating polytope volume", In Proc. ACM Symposium on Computational Geometry 2014, pages 318-325.

### *We develop the R version of VolEsti*  

* We develop a C++ library and a R interface.
* We exclude CGAL and Boost dependecies.
* In folder include we develop the C++ code.
* In folder test you can compile the C++ code using command: g++ vol.cpp -o executable
* In folder R-proj/src we develop the R interface. You can call the C++ library through R using command: Rcpp::sourceCpp('path/to/R-proj/src/vol_R.cpp'). Then function vol_R() is ready to call.
* You have to install Rcpp and RcppEigen libraries.

