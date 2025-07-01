// shakeandbake.cpp
// Usage:
//   ./shakeandbake <cube|simplex|birkhoff|iSDY_1059> <dimension> [epsilon]
// ---------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <limits>

#include <Eigen/Eigen>
#include <boost/random.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/feasible_point.hpp"   
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "random_walks/shake_and_bake_walk.hpp"
#include "generators/known_polytope_generators.h"
#include "diagnostics/scaling_ratio.hpp" 


using NT     = double;
using Kernel = Cartesian<NT>;
using Point  = Kernel::Point;
using RNG    = BoostRandomNumberGenerator<boost::random::mt19937, NT>;
using HPoly  = HPolytope<Point>;
using Walker1 = ShakeAndBakeWalk::Walk<HPoly, RNG>;


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <cube|simplex|birkhoff|iSDY_1059> <dimension> [epsilon]\n";
        return 1;
    }


    std::string shape  = argv[1];
    unsigned    cli_n  = std::stoi(argv[2]);
    NT eps_cli = (argc > 3)
                 ? static_cast<NT>(std::stod(argv[3]))
                 : Walker1::kDefaultEpsilon;      // default value if not manually

    //Generating polytope 
    HPoly P;
    if (shape == "cube")               P = generate_cube<HPoly>(cli_n, false);
    else if (shape == "simplex")       P = generate_simplex<HPoly>(cli_n, false);
    else if (shape == "birkhoff")      P = generate_birkhoff<HPoly>(cli_n);
    else {
        std::cerr << "Unknown polytope type: " << shape << '\n';
        return 1;
    }

    // Walk parameters adjustments
    const unsigned true_dim = P.dimension();
    unsigned walk_len, n_samples, burn_in_iters;

    if (shape == "cube" || shape == "simplex") 
    {
        walk_len      = 50;
        n_samples     = 100000;
        burn_in_iters = 50;
    } 
    else if (shape == "birkhoff") 
    {
        walk_len      = 100 * true_dim;
        n_samples     = 2000 * true_dim;
        burn_in_iters = 10  * true_dim;
    } 
    else {                              
        walk_len      = 20 * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 20 * true_dim;
    }

    std::cout << "Parameters: walk_len="   << walk_len
              << ", n_samples="            << n_samples
              << ", burn_in_iters="        << burn_in_iters
              << " (dim="                  << true_dim
              << ") eps="                  << eps_cli << '\n';

    // Initializing the walk
    RNG rng(true_dim);

    auto [boundary_pt, facet_idx] = compute_boundary_point<Point>(P, rng, eps_cli);
    Walker1 walk1(P, boundary_pt, facet_idx, rng,eps_cli);
    const NT tol = walk1.get_epsilon();                     

    const std::string base = "sb_" + shape + "_" + std::to_string(cli_n);
    std::ofstream out(base + "_run.txt");


    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples1(true_dim, n_samples);
    std::vector<int> facet_id1(n_samples, -1);

    //Burn in 
    for (int i = 0; i < burn_in_iters; ++i)
        walk1.apply(walk_len, rng);

    // Sampling
    for (int i = 0; i < n_samples; ++i) {
        walk1.apply( walk_len, rng);
        const Point& q = walk1.getCurrentPoint();
        samples1.col(i) = q.getCoefficients();

        // File inscription
        for (unsigned d = 0; d < true_dim; ++d)
            out << q[d] << (d + 1 < true_dim ? ' ' : '\n');
    }
    out.close();
    
    std::cout << "Generated " << n_samples << " samples in "
            << walk_len   << " steps each.\n";


    //Scaling ratio test 

    auto [scales, coverage] = scaling_ratio_boundary_test(P, samples1,tol);

    std::cout << "Scaling faktori (x-osa):\n";
    for (double s : scales) {
        std::cout << s << " ";
    }
    std::cout << "\n\nCoverage matrica (svaki red = jedna faseta):\n";
    for (int f = 0; f < coverage.rows(); ++f) {
        std::cout << "Faseta " << f << ": ";
        for (int k = 0; k < coverage.cols(); ++k) {
            double cov = coverage(f, k);
            if (std::isnan(cov))
                std::cout << "NaN ";
            else
                std::cout << cov << " ";
        }
        std::cout << "\n";
    }

    return 0;
}