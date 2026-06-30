// VolEsti (volume computation and sampling library)
// Licensed under GNU LGPL.3, see LICENCE file

// Example: Count linear extensions of a partial order using OrderPolytope
//          volume estimation.
//
// This example reads a poset from a file, constructs the OrderPolytope,
// estimates its volume, and computes #LE = n! * Vol(OrderPolytope).
//
// Usage:
//   ./count_le_orderpolytope <INSTANCE_FILE> <VOLUME_METHOD> [WALK_LENGTH]
//
// Arguments:
//   INSTANCE_FILE  - Poset file (first line: n, then pairs "i j" for a_i <= a_j)
//   VOLUME_METHOD  - "sob", "cg", or "cb"
//   WALK_LENGTH    - Optional (default: auto-computed from dimension)
//
// Examples:
//   ./count_le_orderpolytope instances/chain_4.txt cb
//   ./count_le_orderpolytope instances/antichain_4.txt sob 20
//   ./count_le_orderpolytope instances/diamond_4.txt cg

#include <iostream>
#include <fstream>
#include <chrono>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/hpolytope.h"
#include "misc/poset.h"
#include "misc/misc.h"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/count_linear_extensions.hpp"

#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"


typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
typedef OrderPolytope<Point> ORDERPOLYTOPE;
typedef HPolytope<Point> HPOLYTOPE;
typedef typename ORDERPOLYTOPE::MT MT;
typedef typename ORDERPOLYTOPE::VT VT;


// Count linear extensions using OrderPolytope (optimized)
NT count_le_order_polytope(Poset const& poset, std::string const& algo,
                           unsigned int walk_length)
{
    ORDERPOLYTOPE OP(poset);
    unsigned int d = OP.dimension();
    NT e = 0.1;

    RNGType rng(d);
    NT volume;

    auto start = std::chrono::high_resolution_clock::now();

    if (algo == "sob") {
        volume = volume_sequence_of_balls<CDHRWalk, RNGType>(OP, rng, e, walk_length);
    } else if (algo == "cg") {
        volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(OP, rng, e, walk_length);
    } else {
        volume = volume_cooling_balls<CDHRWalk, RNGType>(OP, rng, e, 2*walk_length).second;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    NT le_count = volume;
    for (NT i = (NT)d; i > 1; i -= 1)
        le_count *= i;

    std::cout << "=== OrderPolytope Method ===" << std::endl;
    std::cout << "  Dimension:            " << d << std::endl;
    std::cout << "  Num hyperplanes:      " << OP.num_of_hyperplanes() << std::endl;
    std::cout << "  Volume:               " << volume << std::endl;
    std::cout << "  Linear extensions:    " << le_count << std::endl;
    std::cout << "  Time (seconds):       " << elapsed << std::endl;

    return le_count;
}


// Count linear extensions using HPolytope (for comparison)
NT count_le_h_polytope(Poset const& poset, std::string const& algo,
                       unsigned int walk_length)
{
    ORDERPOLYTOPE OP(poset);
    unsigned int d = OP.dimension();
    NT e = 0.1;

    // Convert OrderPolytope to HPolytope
    MT A = OP.get_dense_mat();
    VT b_vec = OP.get_vec();
    HPOLYTOPE HP(d, A, b_vec);

    RNGType rng(d);
    NT volume;

    auto start = std::chrono::high_resolution_clock::now();

    if (algo == "sob") {
        volume = volume_sequence_of_balls<CDHRWalk, RNGType>(HP, rng, e, walk_length);
    } else if (algo == "cg") {
        volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, rng, e, walk_length);
    } else {
        volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, rng, e, 2*walk_length).second;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    NT le_count = volume;
    for (NT i = (NT)d; i > 1; i -= 1)
        le_count *= i;

    std::cout << "=== HPolytope Method (comparison) ===" << std::endl;
    std::cout << "  Dimension:            " << d << std::endl;
    std::cout << "  Num hyperplanes:      " << HP.num_of_hyperplanes() << std::endl;
    std::cout << "  Volume:               " << volume << std::endl;
    std::cout << "  Linear extensions:    " << le_count << std::endl;
    std::cout << "  Time (seconds):       " << elapsed << std::endl;

    return le_count;
}


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: ./count_le_orderpolytope INSTANCE_FILE VOLUME_METHOD [WALK_LENGTH]" << std::endl;
        std::cerr << "  INSTANCE_FILE:  poset file (first line: n, then pairs 'i j')" << std::endl;
        std::cerr << "  VOLUME_METHOD:  sob, cg, or cb" << std::endl;
        std::cerr << "  WALK_LENGTH:    optional (default: 10 + d/10)" << std::endl;
        return 1;
    }

    // Parse arguments
    std::string filename(argv[1]);
    std::string algo(argv[2]);

    if (algo != "sob" && algo != "cg" && algo != "cb") {
        std::cerr << "Invalid volume method: " << algo << ". Use sob, cg, or cb." << std::endl;
        return 1;
    }

    // Read poset from file
    std::ifstream data_file(filename);
    if (!data_file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return 1;
    }
    Poset poset = read_poset_from_file(data_file);
    data_file.close();

    std::cout << "Poset loaded: " << poset.num_elem() << " elements, "
              << poset.num_relations() << " relations" << std::endl;

    unsigned int d = poset.num_elem();
    unsigned int walk_length = (argc >= 4) ? std::atoi(argv[3]) : (10 + d / 10);

    std::cout << "Walk length: " << walk_length << std::endl;
    std::cout << "Algorithm:   " << algo << std::endl;
    std::cout << std::endl;

    // Count using OrderPolytope
    NT le_op = count_le_order_polytope(poset, algo, walk_length);

    std::cout << std::endl;

    // Count using HPolytope (for comparison)
    NT le_hp = count_le_h_polytope(poset, algo, walk_length);

    std::cout << std::endl;
    std::cout << "=== Summary ===" << std::endl;
    std::cout << "OrderPolytope LE: " << le_op << std::endl;
    std::cout << "HPolytope LE:     " << le_hp << std::endl;

    return 0;
}
