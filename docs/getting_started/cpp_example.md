# How to create your first example using the C++ interface

Here we give step-by-step instructions for how to estimate the volume of a 3-dimensional cube using `volesti` library.

Write the following C++ code and save it in `volume_example.cpp`

```c++

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_balls.hpp"

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, double, 3> RNGType;
typedef HPolytope<Point> HPolytopeType;

int main() {
	// Generating a 3-dimensional cube centered at origin
	HPolytopeType HP = generate_cube<HPolytopeType>(3, false);
	std::cout<<"Polytope: \n";
	HP.print();
	std::cout<<"\n";

	// Setup parameters for calculating volume
	int walk_len = 10 + HP.dimension()/10;
	double e = 0.1;

	// Calculating volume of the passed polytope
	double volume = volume_cooling_balls
        <BallWalk, RNGType, HPolytopeType>(HP, e, walk_len).second;

    std::cout << "Volume of the cube: " << volume << std::endl;

	return 0;
}
```

Then create a `CMakeList.txt` file with the following text:

```cmake
project( VolEsti-cpp-example )

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

add_definitions(-DDISABLE_NLP_ORACLES)

option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

# Find lpsolve library
find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")
else ()
  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../../external)
  include_directories (BEFORE ../../include)

  # for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard

  add_executable (volume_example volume_example.cpp)
  TARGET_LINK_LIBRARIES(volume_example ${LP_SOLVE})

endif()
```

We will use `cmake` to build the makefile and compile our example.

```bash
cmake .
make
./volume_example
```

The last command will give the following output:

```
Polytope:
 6 3 double
1 0 0 <= 1
0 1 0 <= 1
0 0 1 <= 1
-1 0 0 <= 1
0 -1 0 <= 1
0 0 -1 <= 1

Volume of the cube: 8.22251
```

That is, a 3-dimensional cube defined by 6 inequlity constraints. Our code computes an approximation of the volume.
