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
	double volume = volume_cooling_balls<BallWalk, RNGType, HPolytopeType>(HP, e, walk_len).second;

    std::cout << "Volume of the cube: " << volume << std::endl;

	return 0;
}