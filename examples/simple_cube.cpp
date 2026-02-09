#include "Eigen/Eigen"
#include "generators/known_polytope_generators.h"
#include "volume/volume_sequence_of_balls.hpp"
#include <cmath>
#include <iostream>

int main() {
  typedef double NT;
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef HPolytope<Point> Hpolytope;

  int d = 10;
  std::cout << "Computing volume of a " << d << "-dimensional hypercube..."
            << std::endl;

  // Create a hypercube H-polytope [-1, 1]^d
  Hpolytope HPoly = generate_cube<Hpolytope>(d, false);

  // Compute volume using Sequence of Balls (SOB)
  clock_t start = clock();
  double vol = volume_sequence_of_balls<>(HPoly);
  clock_t end = clock();

  double elapsed = double(end - start) / CLOCKS_PER_SEC;
  double exact_vol = std::pow(2.0, d);
  double error = std::abs(vol - exact_vol) / exact_vol * 100;

  std::cout << "Estimated Volume: " << vol << std::endl;
  std::cout << "Exact Volume:     " << exact_vol << std::endl;
  std::cout << "Relative Error:   " << error << "%" << std::endl;
  std::cout << "Time:             " << elapsed << "s" << std::endl;

  return 0;
}
