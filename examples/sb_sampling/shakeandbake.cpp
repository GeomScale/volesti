// shakeandbake.cpp
// Usage:
//   ./shakeandbake <cube|simplex|birkhoff|iSDY_1059> <dimension> [epsilon]
// ---------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>

#include <Eigen/Eigen>
#include <boost/random.hpp>
#include "cnpy.h"

#include "diagnostics/effective_sample_size.hpp"
#include "diagnostics/interval_psrf.hpp"
#include "diagnostics/print_diagnostics.hpp"

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "random_walks/shake_and_bake_walk.hpp"
#include "generators/known_polytope_generators.h"
#include "misc/print_table.hpp"

using NT     = double;
using Kernel = Cartesian<NT>;
using Point  = Kernel::Point;
using RNG    = BoostRandomNumberGenerator<boost::random::mt19937, NT>;
using HPoly  = HPolytope<Point>;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage : " << argv[0]
                  << " <cube|simplex|birkhoff|iSDY_1059> <dimension> [epsilon]\n";
        return 1;
    }

    std::string shape = argv[1];
    unsigned    cli_n = std::stoi(argv[2]);

    HPoly P;
    if (shape == "cube") {
        P = generate_cube<HPoly>(cli_n, false);
    }
    else if (shape == "simplex") {
        P = generate_simplex<HPoly>(cli_n, false);
    }
    else if (shape == "birkhoff") {
        P = generate_birkhoff<HPoly>(cli_n);
    }
    else if (shape == "iSDY_1059") {
        auto arrA = cnpy::npy_load("A_iSDY_1059.npy");
        auto arrb = cnpy::npy_load("b_iSDY_1059.npy");
        size_t m = arrA.shape[0], n = arrA.shape[1];
        Eigen::MatrixXd A(m,n);
        Eigen::VectorXd b(m);
        double* dataA = arrA.data<double>();
        double* datab = arrb.data<double>();
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j)
                A(i,j) = dataA[i*n + j];
            b(i) = datab[i];
        }
        P = HPoly(n, A, b);
    }
    else {
        std::cerr << "Unknown polytope type: " << shape << "\n";
        return 1;
    }

    const unsigned true_dim = P.dimension();
    unsigned walk_len, n_samples, burn_in_iters;

    //Adaptive tuning - polytope type
    if (shape == "cube") {
        walk_len      = 10  * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 5  * true_dim;
    }
    else if (shape == "simplex") {
        walk_len      = 10  * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 5  * true_dim;
    }
    else if (shape == "birkhoff") {
        walk_len      = 20  * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 5  * true_dim;
    }
    else {
        walk_len      = 20  * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 20  * true_dim;
    }

    std::cout << "Parameters: walk_len=" << walk_len
              << ", n_samples=" << n_samples
              << ", burn_in_iters=" << burn_in_iters
              << " (dim=" << true_dim << ")\n";

    RNG rng(true_dim);
    std::string base  = "sb_" + shape + "_" + std::to_string(cli_n);
    std::string fname = base + "_run.txt";
    std::ofstream out(fname);

    using Walker = ShakeAndBakeWalk::Walk<HPoly, RNG>;
    Walker walk(P, rng, ShakeAndBakeWalk::Running);

    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples(true_dim, n_samples);

    // Burn-in phase for burn_in_iters
    Point tmp0(true_dim);
    for (unsigned i = 0; i < burn_in_iters; ++i) {
        walk.apply(P, tmp0, walk_len, rng);
    }

    // Actual sampling
    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < n_samples; ++i) {
        Point tmp(true_dim);
        walk.apply(P, tmp, walk_len, rng);
        const Point& q = walk.getCurrentPoint();
        for (unsigned j = 0; j < true_dim; ++j) {
            out << q[j] << (j + 1 < true_dim ? ' ' : '\n');
            samples(j,i) = q[j];
        }
    }
    out.close();
    auto t1 = std::chrono::high_resolution_clock::now();

    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Mode Shake-and-Bake Running: "
              << n_samples << " samples (walk_len=" << walk_len
              << ", burn_in=" << burn_in_iters << ") generated in "
              << secs << " s\n";

    //Diagnostics
    unsigned min_ess = 0;
    std::cout << "== Diagnostics per dimension ==\n";
    print_diagnostics<NT, Eigen::VectorXd, decltype(samples), std::ostream>
        (samples, min_ess, std::cout);
    std::cout << "Minimum ESS (for thinning): " << min_ess << "\n\n";

    // Full polytope "uniformity test"
    Eigen::MatrixXd A = P.get_mat();  
    Eigen::VectorXd b = P.get_vec();  
    Eigen::VectorXd centroid = samples.rowwise().mean();  

    std::string tfile = base + "_tvals.txt";
    std::ofstream out_t(tfile);

    
    Eigen::VectorXd rhs = b - A * centroid;  

    for (unsigned i = 0; i < n_samples; ++i) {

        Eigen::VectorXd X = samples.col(i);     
        Eigen::VectorXd diff = X - centroid;    

        Eigen::VectorXd denom = A * diff;       

        double t_i = 1.0;
        for (int k = 0; k < denom.size(); ++k) {
            if (denom[k] > 0) {
                double tlim = rhs[k] / denom[k];
                if (tlim < t_i) t_i = tlim;
            }
        }
        if (t_i < 0.0) t_i = 0.0;
        if (t_i > 1.0) t_i = 1.0;

        out_t << t_i << "\n";
    }
    out_t.close();
    std::cout << "Wrote to " << tfile << "\n";

    return 0;
}