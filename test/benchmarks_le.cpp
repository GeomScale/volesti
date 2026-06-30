// VolEsti (volume computation and sampling library)
// Licensed under GNU LGPL.3, see LICENCE file

// Benchmarks for linear extension counting:
//   - OrderPolytope vs HPolytope based volume estimation
//   - Different volume algorithms (SOB, CG, CB)
//   - Varying poset sizes and densities

#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <random>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/hpolytope.h"
#include "misc/poset.h"

#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"


typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
typedef OrderPolytope<Point> ORDERPOLYTOPE;
typedef HPolytope<Point> HPOLYTOPE;
typedef typename ORDERPOLYTOPE::MT MT;
typedef typename ORDERPOLYTOPE::VT VT;


// Generate a random poset with n elements and approximately density * n*(n-1)/2 relations
Poset generate_random_poset(unsigned int n, double density, unsigned int seed = 42) {
    typedef typename Poset::RT RT;
    typedef typename Poset::RV RV;

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Generate a random DAG by only adding edges from lower to higher indexed nodes
    RV relations;
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = i + 1; j < n; ++j) {
            if (dis(gen) < density) {
                relations.push_back(RT(i, j));
            }
        }
    }

    return Poset(n, relations);
}


struct BenchmarkResult {
    std::string method;
    std::string algo;
    unsigned int n;
    unsigned int num_relations;
    NT volume;
    NT le_count;
    double time_seconds;
};


// Run a single benchmark
template <typename Polytope, typename WalkPolicy>
BenchmarkResult run_benchmark(Polytope& P, unsigned int n, unsigned int num_relations,
                               std::string const& method, std::string const& algo,
                               NT error, unsigned int walk_length) {
    RNGType rng(P.dimension());
    BenchmarkResult result;
    result.method = method;
    result.algo = algo;
    result.n = n;
    result.num_relations = num_relations;

    auto start = std::chrono::high_resolution_clock::now();

    if (algo == "SOB") {
        result.volume = volume_sequence_of_balls<WalkPolicy, RNGType>(P, rng, error, walk_length);
    } else if (algo == "CG") {
        result.volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(P, rng, error, walk_length);
    } else {
        result.volume = volume_cooling_balls<WalkPolicy, RNGType>(P, rng, error, walk_length).second;
    }

    auto end = std::chrono::high_resolution_clock::now();
    result.time_seconds = std::chrono::duration<double>(end - start).count();

    result.le_count = result.volume;
    for (NT i = (NT)n; i > 1; i -= 1)
        result.le_count *= i;

    return result;
}


void print_header() {
    std::cout << std::left
              << std::setw(15) << "Method"
              << std::setw(8) << "Algo"
              << std::setw(6) << "n"
              << std::setw(10) << "Relations"
              << std::setw(16) << "Volume"
              << std::setw(16) << "LE Count"
              << std::setw(12) << "Time(s)"
              << std::endl;
    std::cout << std::string(83, '-') << std::endl;
}


void print_result(BenchmarkResult const& r) {
    std::cout << std::left
              << std::setw(15) << r.method
              << std::setw(8) << r.algo
              << std::setw(6) << r.n
              << std::setw(10) << r.num_relations
              << std::setw(16) << std::scientific << std::setprecision(4) << r.volume
              << std::setw(16) << r.le_count
              << std::setw(12) << std::fixed << std::setprecision(3) << r.time_seconds
              << std::endl;
}


int main(int argc, char* argv[])
{
    std::cout << "======================================================================" << std::endl;
    std::cout << "   Linear Extension Counting: OrderPolytope vs HPolytope Benchmarks   " << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    NT error = 0.5;  // relatively loose for benchmarks

    // Test configurations: (n, density)
    std::vector<std::pair<unsigned int, double>> configs = {
        {4, 0.0},   // antichain: 4! = 24 LEs
        {4, 0.3},   // sparse poset
        {4, 1.0},   // chain: 1 LE
        {6, 0.0},   // antichain: 6! = 720 LEs
        {6, 0.3},   // sparse poset
        {8, 0.2},   // medium sparse
        {10, 0.15}, // larger sparse
    };

    std::vector<std::string> algos = {"SOB", "CB"};

    print_header();

    for (auto const& config : configs) {
        unsigned int n = config.first;
        double density = config.second;

        Poset poset = generate_random_poset(n, density);
        unsigned int walk_length = 10 + n / 10;

        ORDERPOLYTOPE OP(poset);
        MT A = OP.get_dense_mat();
        VT b_vec = OP.get_vec();
        HPOLYTOPE HP(n, A, b_vec);

        for (auto const& algo : algos) {
            // OrderPolytope
            ORDERPOLYTOPE OP2(poset);
            auto r1 = run_benchmark<ORDERPOLYTOPE, CDHRWalk>(
                OP2, n, poset.num_relations(), "OrderPolytope", algo, error, walk_length);
            print_result(r1);

            // HPolytope
            HPOLYTOPE HP2(n, A, b_vec);
            auto r2 = run_benchmark<HPOLYTOPE, CDHRWalk>(
                HP2, n, poset.num_relations(), "HPolytope", algo, error, walk_length);
            print_result(r2);
        }

        std::cout << std::string(83, '-') << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Benchmark complete." << std::endl;

    return 0;
}
