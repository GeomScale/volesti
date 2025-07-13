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

using Walker1 = BilliardShakeAndBakeWalk::Walk<HPoly, RNG>;

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
                 : Walker1::kDefaultEps;      

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
              << ") eps="                << eps_cli << '\n';

    RNG rng(true_dim);
    auto [boundary_pt, facet_idx] = compute_boundary_point<Point>(P, rng, eps_cli);

    Walker1 walk1(P, boundary_pt, facet_idx, rng, eps_cli, 3.5);
    const NT tol = walk1.get_epsilon();

    const std::string base = "billiard_sb_" + shape + "_" + std::to_string(cli_n);
    std::ofstream out(base + ".txt");

    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples1(true_dim, n_samples);

    for (int i = 0; i < burn_in_iters; ++i)
        walk1.apply(walk_len, rng);

    for (int i = 0; i < n_samples; ++i) {
        walk1.apply(walk_len, rng);
        const Point& q = walk1.getCurrentPoint();
        samples1.col(i) = q.getCoefficients();

        for (unsigned d = 0; d < true_dim; ++d)
            out << q[d] << (d + 1 < true_dim ? ' ' : '\n');
    }
    out.close();

    std::cout << "Generated " << n_samples << " samples in "
              << walk_len   << " steps each.\n";

    return 0;
}
