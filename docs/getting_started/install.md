Installation Guide
==================

## C++ Interface

### Compile C++ sources and run tests
---

To compile the C++ code you need the [lp_solve](http://lpsolve.sourceforge.net/5.5/) library. For example, for Unix/Linux you need `liblpsolve55.so`. This is available from the library's [webpage](http://lpsolve.sourceforge.net/5.5/) as well as a package in several linux distributions e.g. [debian](https://packages.debian.org/stretch/liblpsolve55-dev) `sudo apt-get install lp-solve`.

You have to specify the path to `liblpsolve55.so/dll/dylib`, by running, in folder test:

```bash
mkdir -p test/build && cd test/build
cmake -DLP_SOLVE=_PATH_TO_LIB_FILE_ ..
make
```
For example:  `-DLP_SOLVE=/usr/lib/lpsolve/liblpsolve55.so`

You can run the tests by `cmake test` or `ctest -jK` where `K` the number of `CPU` threads. By adding the option `--verbose` to `ctest` you get more information about the tests, *e.g.* time per test, volume computed and the name of the polytope or convex body.

### Development environment from Docker container
---
Optionally, it is possible to setup a docker contaner with development environment. To get started with docker, see
[here](https://docs.docker.com/get-started/). Below, here is how a Dockerfile can be written
```
FROM ubuntu:18.04

RUN apt-get update && apt-get install -y g++ cmake lp-solve && \
    rm -rf /var/lib/apt/lists/*
```

The user should create a file `Dockerfile.dev` with above content inside a temporary folder (e.g. `docker`) then create an image and run a container:

```bash
# build docker image
cd docker
docker build -t volesti:dev -f Dockerfile.dev .
# check built image
docker images | grep volesti
# run a container in an interactive mode from volesti source folder
docker run -it -v $PWD:/volesti -w /volesti --name=volesti-dev volesti:dev /bin/bash
```

## R Interface

### Install Rcpp package
---

1. Install package-dependencies: ``Rcpp``, ``RcppEigen``, ``BH``.

2. Then use ``devtools`` package to install ``volesti`` Rcpp package. From terminal go to folder ``/root/R-proj`` and run in terminal:

```bash
    Rscript -e 'Rcpp::compileAttributes()'
    R CMD INSTALL --no-multiarch --with-keep.source .
```

3. You can use Rstudio as well to open ``volesti.Rproj`` and then click `build source Package` and then `Install and Restart` in `Build` at the menu bar.

### Generate CRAN version
---

To generate the CRAN version of the R package follow the instructions below:

1. From the command line navigate to folder ``/cran_gen``. Then Run:

```r
    source('genCRANpkg.R')
```

2. Open ``genCRANpkg.R`` script with `Rstudio` and run it.

### Run volesti from R
---

* The main function is ``volume()``. It can be used to approximate the volume of a convex polytope given as a set of linear inequalities or a set of vertices (d-dimensional points) or as a Minkowski sum of segments (zonotope). There are three algorithms that can be used (``SequenceOfBalls``, ``CoolingGaussian`` and ``CoolingBalls``).
* The function ``sample_points()`` can be used to sample points from a convex polytope with uniform or spherical gaussian target distribution.
* The function ``round_polytope()`` can be used to round a convex polytope.
* The function ``rand_rotate()`` can be used to apply a random rotation to a convex polytope.

For more details you can read the documentation in folder ``/inst/doc``.


## Python Interface

A ``python`` interface is available from the package [dingo](https://github.com/GeomScale/dingo).


