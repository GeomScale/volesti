#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "known_polytope_generators.h"
#include "random_walks/random_walks.hpp"

#include "vpolytope.h"
#include <iostream>
#include <fstream>
#include "misc.h"
#include <limits>

#include "volume_cooling_gaussians.hpp"
#include "volume_cooling_balls.hpp"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef VPolytope <Point> VPOLYTOPE;

void calculateVolumes(VPOLYTOPE &VP) {
	// Setup parameters for calculating volume
	int walk_len = 10 + VP.dimension()/10;
	NT e=0.1;

	// Calculating volume of the passed polytope
	NT volume2 = volume_cooling_gaussians<GaussianBallWalk, RNGType>(VP, e, walk_len);

	std::cout<<"\t Using Cooling Gaussians method: "<<volume2<<"\n";
}


int main(int argc, char* argv[]) {

	// Generating a 4-dimensional VPolytope
	VPOLYTOPE VP1 = generate_cross<VPOLYTOPE>(4, true);
	std::cout<<"Polytope VP1: \n";
	VP1.print();
	std::cout<<"\n";

	std::cout<<"Volume of VP1: \n";
	calculateVolumes(VP1);

	return 0;
}