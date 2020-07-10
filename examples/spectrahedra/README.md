# Examples for Spectrahedra

## Table of contents
1. [Compilation](#compilation)
    1. [Dependencies](#dependencies)
2. [Examples](#examples)
    1. [Read/write SDPA format files - read_write_sdpa_file.cpp](#readwrite-sdpa-format-files---read_write_sdpa_filecpp)
    2. [Sample with HMC, Boltzmann distribution - boltzmann_hmc_walk.cpp](#sample-with-hmc-boltzmann-distribution---boltzmann_hmc_walkcpp)
    3. [Randomized SDP Solver - solve_sdp.cpp](#randomized-sdp-solver---solve_sdpcpp)


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

And copy the folder external back in folder examples/spectrahedra:

<br>

## Examples
### Read/write SDPA format files - read_write_sdpa_file.cpp

In this example, we will read a semidefinite program from a SDPA format input file, print it
and then write it to a new SDPA format file. Run the example with:

```bash
./read_write_sdpa_file
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


### Sample with HMC, Boltzmann distribution - boltzmann_hmc_walk.cpp

In this example, we will sample a spectrahedron under the Boltsmann distribution e^(-c*x/T), using
the hamiltonian monte carlo random walk with reflections. We will read the spectrahedron as
in [readWriteSdpaFile.cpp](#readwrite-sdpa-format-files---readwritesdpafilecpp). Run the example with:

```bash
./boltzmann_hmc_walk
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
HmcWalk hmcWalk(settings);
hmcWalk.apply(spectrahedron, initialPoint, pointsNum, points);
```

- spectrahedron: instance of class Spectrahedron
- initialPoint: an interior point in the spectrahedron
- pointsNum: how many points to sample
- points: a list to return the samples


### Randomized SDP Solver - solve_sdp.cpp

In this example, we will solve a semidefinite program. We will read the program 
as in [read_write_sdpa_file.cpp](#readwrite-sdpa-format-files---read_write_sdpa_filecpp). Run the example with:

```bash
./solve_sdp
```

#### Code Explanation
To use the solver, first we declare some parameters:

```bash
SimulatedAnnealingSettings<Point> settings(rel_error);
```

Actually, we can further customize the algorithm. The full settings definition is:

```bash
SimulatedAnnealingSettings<Point> settings(rel_error, walkLength, maxNumSteps, k)
```

- rel_error: The desired relative error.
- walkLength: Default and recommended is 1. This solver uses the [HMC random walk](#sample-with-hmc-boltzmann-distribution---boltzmann_hmc_walkcpp).
  How many points the walk will "burn" before returning a sample.
- maxNumSteps: Default is -1 (infinite). How many steps we will allow the algorithm.
- k: Default is 0.5. Lower values may achieve faster convergence.

Next we can solve the program:

```bash
NT min = solve_sdp(spectrahedron, objFunction, settings, initialPoint, sol ,verbose);
```

- spectrahedron: Instance of class Spectrahedron (a linear matrix inequality).
- objFunction: The objective function of the program.
- Settings: As above.
- initialPoint: An interior point in the spectrahedron.
- min: The estimated minimum value
- sol: At which point in the spectrahedron (returned by the solver)
- verbose: If true, print useful information.
