#include "Eigen/Eigen"
#include "generators/known_polytope_generators.h"
#include "volume/volume_sequence_of_balls.hpp"
#include <iostream>

// Mission 1: The Simple Cube (Fixed)
// Computes the volume of a 10-dimensional hypercube.

int main() {
  // 1. Define Types (Standard Volesti Setup)
  typedef double NT;
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef HPolytope<Point> Hpolytope;

  // 2. Define Dimension
  int d = 10;
  std::cout << "Computing volume of a " << d << "-dimensional hypercube..."
            << std::endl;

  // 3. Create Hypercube
  // generate_cube<Type>(dimension, is_V_polytope)
  Hpolytope HPoly = generate_cube<Hpolytope>(d, false);

  // 4. Compute Volume
  // Using Sequence of Balls (SOB) algorithm, which is robust.
  // The exact volume is 2^10 = 1024.

  // We clock the computation for fun
  double tstart = (double)clock() / (double)CLOCKS_PER_SEC;

  double vol = volume_sequence_of_balls<>(HPoly);

  double time = (double)clock() / (double)CLOCKS_PER_SEC - tstart;

  // 5. Output Results
  std::cout << "Estimated Volume: " << vol << std::endl;
  std::cout << "Exact Volume:     " << std::pow(2.0, d) << std::endl;
  std::cout << "Time Taken:       " << time << " seconds" << std::endl;

  double error = std::abs(vol - std::pow(2.0, d)) / std::pow(2.0, d) * 100;
  std::cout << "Error:            " << error << "%" << std::endl;

  return 0;
}
