###  R Interface
------------

####  Install Rcpp package  
 
* Install package-dependencies: `Rcpp`, `RcppEigen`, `BH`, `lpSolveAPI`.  

1. Then use devtools package to install `volesti` Rcpp package. In folder /root/R-proj Run:
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
