#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "orderpolytope.h"
 #include "random_walks/gaussian_accelerated_billiard_walk.hpp"
 #include "generators/boost_random_number_generator.hpp"
#include "volume_cooling_ellipsoids.hpp"

#include <iostream>
#include <fstream>
#include "misc.h"
#include "poset.h"


typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef OrderPolytope<Point> ORDER_POLYTOPE;


NT calculateLinearExtension(ORDER_POLYTOPE const& OP) {
    // Setup parameters for calculating volume and rounding
    unsigned int d = OP.dimension();
    unsigned int walk_len = 10 + d/10;
    NT e=0.1;

    NT volume = volume_cooling_ellipsoids<GaussianAcceleratedBilliardWalk, RNGType>(OP, e, 2*walk_len).second;

    // multiplying by d factorial, d = number of elements
    for(NT i=(NT)d; i>1; i-=1) {
        volume = volume * i;
    }

    return volume;
}


/**

 Usage: ./volesti_lecount_order_polytope INSTANCE

 example:
    ./volesti_lecount_order_polytope instances/bipartite_0.5_008_0.txt

*/
int main(int argc, char* argv[]) {
    std::string filename (argv[1]);
    std::pair<bool, Poset> read_poset = read_poset_from_file_adj_matrix(filename);

    if( !read_poset.first ) {
        std::cerr << "error in reading data from file" << std::endl;
        return -1;
    }

    ORDER_POLYTOPE OP(read_poset.second);
    OP.normalize();
    std::cout << calculateLinearExtension(OP) << std::endl;
    return 0;
}
