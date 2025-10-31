// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Eigen>
#include <boost/random.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/feasible_point.hpp"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "random_walks/billiard_shake_and_bake_walk.hpp"
#include "generators/known_polytope_generators.h"
#include "diagnostics/scaling_ratio.hpp"

using NT     = double;
using Kernel = Cartesian<NT>;
using Point  = Kernel::Point;
using RNG    = BoostRandomNumberGenerator<boost::random::mt19937, NT>;
using HPoly  = HPolytope<Point>;
using RMode  = BilliardShakeAndBakeWalk::ReflectionMode;

using Walker1 = BilliardShakeAndBakeWalk::Walk<HPoly, RNG>;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <cube|simplex|birkhoff> <dimension> [epsilon]\n";
        return 1;
    }

    std::string shape  = argv[1];
    unsigned    cli_n  = std::stoi(argv[2]);

    // Defaults
    int nr_cli = -1; // 0 => auto = ceil(sqrt(dim))

    // argv[3] => nr (ako postoji)
    if (argc > 3) nr_cli = std::stoi(argv[3]);

    NT eps_cli = (argc > 4)
                 ? static_cast<NT>(std::stod(argv[4]))
                 : Walker1::kDefaultEpsilon; 

    // POSLE
    RMode rmode = RMode::InverseExponential;
    if (argc > 5) {
        std::string dist = argv[5];
        if (dist == "uniform")      rmode = RMode::Uniform;
        else if (dist == "inverseexp") rmode = RMode::InverseExponential;
        else {
            std::cerr << "Unknown mode: " << dist << " (use: uniform | inverseexp)\n";
            return 1;
        }
    }


    HPoly P;
    if (shape == "cube")      P = generate_cube<HPoly>(cli_n, false);
    else if (shape == "simplex") P = generate_simplex<HPoly>(cli_n, false);
    else if (shape == "birkhoff") P = generate_birkhoff<HPoly>(cli_n);
    else {
        std::cerr << "Unknown polytope type: " << shape << '\n';
        return 1;
    }

    const unsigned true_dim   = P.dimension();
    unsigned       walk_len, n_samples, burn_in_iters;

    int mode = (shape == "cube" || shape == "simplex") ? 0
             : (shape == "birkhoff")             ? 1
             : 2;

    switch (mode) {
        case 0:
            walk_len      = 20  * true_dim;
            n_samples     = 500 * true_dim;
            burn_in_iters = 5   * true_dim;
            break;
        case 1:
            walk_len      = 100 * true_dim;
            n_samples     = 2000 * true_dim;
            burn_in_iters = 10  * true_dim;
            break;
        default:
            walk_len      = 20  * true_dim;
            n_samples     = 100 * true_dim;
            burn_in_iters = 20  * true_dim;
            break;
    }

    std::cout << "Parameters: walk_len="   << walk_len
              << ", n_samples="          << n_samples
              << ", burn_in_iters="      << burn_in_iters
              << " (dim="                << true_dim
              << ", nr="                 << nr_cli
              << ") eps="                << eps_cli << '\n';

    RNG rng(true_dim);
    auto [boundary_pt, facet_idx] = compute_boundary_point<Point>(P, rng, eps_cli);

    Walker1 walk1(P, boundary_pt, rng, facet_idx, nr_cli, eps_cli, rmode);
    const NT tol = walk1.get_epsilon();

    const std::string base = "billiard_sb_" + shape + "_" + std::to_string(cli_n)+"_" + std::to_string(nr_cli);
    std::ofstream out(base + ".txt");

    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples1(true_dim, n_samples);

    //Burn in 
    for (int i = 0; i < burn_in_iters; ++i)
        walk1.apply(P, walk_len, rng);

    // Sampling
    for (int i = 0; i < n_samples; ++i) {
        walk1.apply(P, walk_len, rng);
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
    auto [scales, coverage, max_dev, avg_dev] = scaling_ratio_boundary_test(P, samples1,tol);

    std::cout << "Scaling factors:\n";
    for (double s : scales) {
        std::cout << s << " ";
    }
    std::cout << "\n\nCoverage matrix (each row = one facet):\n";
    for (int f = 0; f < coverage.rows(); ++f) {
        std::cout << "Facet " << f << ": ";
        for (int k = 0; k < coverage.cols(); ++k) {
            double cov = coverage(f, k);
            if (std::isnan(cov))
                std::cout << "NaN ";
            else
                std::cout << cov << " ";
        }
        std::cout << "\n";
    }

    //Uniformity deviation analysis 
    std::cout << "\n";
    std::cout << "Facet        Max deviation (%)        Avg deviation (%)\n";
    for (int f = 0; f < max_dev.size(); ++f) {
        std::cout << std::setw(6) << f << " "
                  << std::fixed << std::setprecision(2)
                  << std::setw(18) << max_dev[f] << " "
                  << std::setw(22) << avg_dev[f]
                  << "\n";
    }
    return 0;
}