# Volume computation and sampling

## About
The `volesti` package provides [R](https://www.r-project.org/) with functions for volume estimation and sampling. In particular, it provides an R interface for the C++ library [**volesti**](https://github.com/GeomScale/volesti).

`volesti` computes approximations of volume of polytopes given as a set of points or linear inequalities or as a Minkowski sum of segments (zonotopes). There are algorithms for volume approximation as well as algorithms for sampling, rounding and rotating polytopes. Last but not least, `volesti` provides implementations of geometric algorithms to compute the score of a portfolio given asset returns and to detect financial crises in stock markets.

##  Download and install

* The latest stable version is available from CRAN.
* The latest development version is available on Github `www.github.com/GeomScale/volesti`

* Install `volesti` by running:
```
install.packages("volesti")
```
* The package-dependencies are: `Rcpp`, `RcppEigen`, `BH`.

## Documentation

* [Using the R Interface](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/r_interface.md)
* [Wikipage with Tutorials and Demos](https://github.com/GeomScale/volesti/wiki)
* [Tutorial given to PyData meetup](https://vissarion.github.io/tutorials/volesti_tutorial_pydata.html)


## Credits

* [Contributors and Package History](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/credits.md)
* [List of Publications](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/publications.md)

Copyright (c) 2012-2020 Vissarion Fisikopoulos
Copyright (c) 2018-2020 Apostolos Chalkis

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.

