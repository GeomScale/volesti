# Examples for Spectrahedra

## Table of contents
1. [Compilation](#compilation)
    1. [Dependencies](#dependencies)
2. [Examples](#examples)
    1. [Example 1: Read/write SDPA format files](#example-1-readwrite-sdpa-format-files)
    2. [Example 2: Sample with HMC, Boltzmann distribution](#example-2-sample-with-hmc-boltzmann-distribution)

## Compilation
In folder examples, first run cmake, to create the makefile:

```bash
cmake .
```

Then, in folder examples/spectrahedra compile using the makefile:

```bash
make
```

### Dependencies
To compile some programs in this folder, we need the libraries openblas, lapack and arpack. If you want to compile 
using the provided cmakelists file, follow the next steps to install and link them. 


First we will need a [Fortran compiler](https://gcc.gnu.org/wiki/GFortran) for GCC. In linux:
```bash
sudo apt install gfortran
```

You may have to edit the path in following line in examples/spectrahedra/CMakeLists.txt

```bash
FIND_LIBRARY(GFORTRAN_LIB NAMES libgfortran.so PATHS /usr/lib/gcc/x86_64-linux-gnu/8/)
```

to point to your installation.
Then we can install the openblas, lapack and arpack libraries (lapack is included in openblas). 
In the folder "examples", clone this repo:

```bash
git clone https://github.com/m-reuter/arpackpp
cd arpackcpp
```

It has two scripts that should easily install the libraries:

```bash
./install-openblas.sh
./install-arpack-ng.sh
```

And copy the folder external back in folder examples:

```bash
 cp -r external ../
```

<br>

## Examples
### Example 1: Read/write SDPA format files

In this example, we will read a semidefinite program from a SDPA format input file, print it
and then write it to a new SDPA format file. Run the example with:

```bash
./readWriteSdpaFile
```

The input file is data/sdp_n2m3.txt. It contains a semidefinite program in SDPA format. A semidefinite program
(dual form) consists of a linear matrix inequality (describing a spectrahedron) of the form A1 + ... + xn An - A0 >= 0,
where >= 0 denotes positive semidefiniteness. For details in the SDPA format you can search the official manual.

The contents of the file are:

```bash
2
1
3
1 1
-1  0  0
 0 -2  1
 0  1 -2
 1 -0 -0
-0 -0 -1
-0 -1 -0
-0 -0  1
-0 -0 -0
 1 -0 -0
```

 
It represents a spectrahedron in 2 dimensions, described by a linear matrix inequality with
3x3 matrices, and a linear objective function. 

##### Explanation line by line
- 2 : The number of dimensions. The vector of the objective function is of length 2 and we need 3 matrices for 
   the linear matrix inequality, i.e. x A1 + y A2 - A0 >= 0.
- 1 : Please lookup SDPA format for more details. Currently, support is only for value 1.
- 3 : The size of the matrices, i.e. 3x3
- 1 1 : The vector of the objective function, i.e. (1,1) or f(x,y) = x + y
- -1  0  0: The first row of A0
- 0 -2  1: The second row of A0
- 0  1 -2: The third row of A0
- 1 -0 -0: The first row of A1
- and so on, till all 3 matrices are defined


### Example 2: Sample with HMC, Boltzmann distribution

In this example, we will sample a spectrahedron under the Boltsmann distribution e^(-c*x/T), using
the hamiltonian monte carlo random walk with reflections. We will read the spectrahedron as
in [Example 1](#example-1-readwrite-sdpa-format-files). Run the example with:

```bash
./boltzmannHmcWalk
```

#### Code Explanation
In boltzmannHmcWalk.cpp, to use the random walk first we need to declare some parameters:

```bash
HmcWalkSettings settings(walkLength, randomNumberGenerator, objFunction, temperature, diameter);
```

- walkLength: how many points the walk will "burn" before returning a sample
- randomNumberGenerator: a class that generates random numbers
- objFunction: the vector c in the boltzmann distribution e^(-c*x/T)
- temperature: T in e^(-c*x/T)
- diameter: diameter of the spectrahedron; can estimate it with a heuristic - method of class Spectrahedron

and then we can sample the spectrahedron

```bash
hmcWalk.apply(spectrahedron, initialPoint, pointsNum, points);
```

- sepctrahedron: instance of class Spectrahedron
- initialPoint: an interior point in the spectrahedron
- pointsNum: how many points to sample
- points: a list to return the samples

