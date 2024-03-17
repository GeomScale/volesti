
#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "hpolytope.h"
#include "known_polytope_generators.h"

#include "random_walks/random_walks.hpp"

#include "volume_sequence_of_balls.hpp"
#include "volume_cooling_gaussians.hpp"
#include "volume_cooling_balls.hpp"

#include <iostream>
#include <fstream>
#include "misc.h"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef HPolytope <Point> HPOLYTOPE;

void calculateVolumes(HPOLYTOPE &HP) {
	// Setup parameters for calculating volume
	int walk_len = 10 + HP.dimension()/10;
	NT e=0.1;

	// Calculating volume of the passed polytope
	NT volume = volume_cooling_gaussians<GaussianHamiltonianMonteCarloExactWalk, RNGType>(HP, e, walk_len);

	std::cout<<"\t Using Cooling Gaussians method: "<<volume<<"\n";
	
}


int main(int argc, char* argv[]) {

	// Generating a 3-dimensional cube centered at origin
	HPOLYTOPE HP = generate_cube<HPOLYTOPE>(3, false);
	std::cout<<"Polytope HP: \n";
	HP.print();
	std::cout<<"\n";

	std::cout<<"Volume of HP: \n";
	calculateVolumes(HP);

	return 0;
}