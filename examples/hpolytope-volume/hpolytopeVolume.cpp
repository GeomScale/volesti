// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Contributed and/or modified by Vaibhav Thakkar

// Licensed under GNU LGPL.3, see LICENCE file

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

void calculateVolumes(const HPOLYTOPE &HP) {
	// Setup parameters for calculating volume
	int walk_len = 10 + HP.dimension()/10;
	NT e=0.1;

	// Calculating volume of the passed polytope
	NT volume1 = volume_cooling_balls<BallWalk, RNGType, HPOLYTOPE>(HP, e, walk_len);
	NT volume2 = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, e, walk_len);
	NT volume3 = volume_sequence_of_balls<BallWalk, RNGType>(HP, e, walk_len);

	std::cout<<"\t Using Cooling Balls method: "<<volume1<<"\n";
	std::cout<<"\t Using Cooling Gaussians method: "<<volume2<<"\n";
	std::cout<<"\t Using Sequence of Balls method: "<<volume3<<"\n";
}


int main(int argc, char* argv[]) {
	// Generating a 4-dimensional cube centered at origin
	HPOLYTOPE HP1 = gen_cube<HPOLYTOPE>(4, false);
	std::cout<<"Polytope HP1: \n";
	HP1.print();
	std::cout<<"\n";

	std::cout<<"Volume of HP1: \n";
	calculateVolumes(HP1);


	// Reading a polytope from ine file
	std::string fileName("data/cube10.ine");
	std::cout<<"Reading input from file..."<<std::endl;
	std::ifstream inp;
	std::vector<std::vector<NT> > Pin;
	inp.open(fileName, std::ifstream::in);
	read_pointset(inp,Pin);

	HPOLYTOPE HP2;
	HP2.init(Pin);
	std::cout<<"Polytope HP2: \n";
	HP2.print();
	std::cout<<"\n";

	std::cout<<"Volume of HP2: \n";
	calculateVolumes(HP2);

	return 0;
}