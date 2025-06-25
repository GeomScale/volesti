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
#include "cnpy.h"

#include "diagnostics/effective_sample_size.hpp"
#include "diagnostics/interval_psrf.hpp"
#include "diagnostics/print_diagnostics.hpp"

#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/feasible_point.hpp"   
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
using Walker = ShakeAndBakeWalk::Walk<HPoly, RNG>;

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
                 : Walker::kDefaultEpsilon;      // default value if not manually

    //Generating polytope 
    HPoly P;
    if (shape == "cube")               P = generate_cube<HPoly>(cli_n, false);
    else if (shape == "simplex")       P = generate_simplex<HPoly>(cli_n, false);
    else if (shape == "birkhoff")      P = generate_birkhoff<HPoly>(cli_n);
    else if (shape == "iSDY_1059") // this is a random biological polytope from py just ignore :)
    {
        auto arrA = cnpy::npy_load("A_iSDY_1059.npy");
        auto arrb = cnpy::npy_load("b_iSDY_1059.npy");
        const size_t m = arrA.shape[0], n = arrA.shape[1];

        Eigen::MatrixXd A(m, n);
        Eigen::VectorXd b(m);

        double* dataA = arrA.data<double>();
        double* datab = arrb.data<double>();
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j)
                A(i, j) = dataA[i * n + j];
            b(i) = datab[i];
        }
        P = HPoly(n, A, b);
    }
    else {
        std::cerr << "Unknown polytope type: " << shape << '\n';
        return 1;
    }

    // Walk parameters adjustments
    const unsigned true_dim = P.dimension();
    unsigned walk_len, n_samples, burn_in_iters;

    if (shape == "cube" || shape == "simplex") 
    {
        walk_len      = 20;
        n_samples     = 100000;
        burn_in_iters = 5;
    } 
    else if (shape == "birkhoff") 
    {
        walk_len      = 20 * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 5  * true_dim;
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

    Walker walk(P, boundary_pt, facet_idx, rng,eps_cli);
    const NT tol = walk.get_epsilon();                     

    const std::string base = "sb_" + shape + "_" + std::to_string(cli_n);
    std::ofstream out(base + "_run.txt");

    const auto A_full = P.get_mat();
    const auto b_full = P.get_vec();
    const size_t m    = A_full.rows();

    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples(true_dim, n_samples);
    std::vector<int> facet_id(n_samples, -1);

    //Burni in 
    for (int i = 0; i < burn_in_iters; ++i)
        walk.apply(walk_len, rng);

    // Sampling
    for (int i = 0; i < n_samples; ++i) {
        walk.apply( walk_len, rng);
        const Point& q = walk.getCurrentPoint();

        Eigen::VectorXd qv(true_dim);
        for (int d = 0; d < true_dim; ++d) {
            qv[d]  = q[d];
            samples(d, i) = q[d];
        }

        const auto Aq = A_full * qv;
        for (size_t k = 0; k < m; ++k) {
            if (std::abs(Aq[k] - b_full[k]) < tol) {
                facet_id[i] = static_cast<int>(k);
                break;
            }
        }

        // File inscription
        for (unsigned d = 0; d < true_dim; ++d)
            out << q[d] << (d + 1 < true_dim ? ' ' : '\n');
    }
    out.close();

    std::cout << "Generated " << n_samples << " samples in "
              << walk_len   << " steps each.\n";

    // Scaling ratio test
    constexpr double min_ratio = 0.01;
    std::ofstream cov_out(base + "_coverage.txt");
    cov_out << "\nScaling coverage by facet (skip <" << min_ratio << "):\n";

    for (int f = 0; f < static_cast<int>(m); ++f) 
    {

        // Samples S on facet f
        std::vector<int> S;
        for (unsigned i = 0; i < n_samples; ++i)
            if (facet_id[i] == f) S.push_back(i);

        const double ratio = static_cast<double>(S.size()) / n_samples;
        if (ratio < min_ratio) {
            cov_out << "facet " << f << " skipped (" << ratio << ")\n";
            continue;
        }

        // Finding the center
        Eigen::VectorXd p = Eigen::VectorXd::Zero(true_dim);
        for (int idx : S) p += samples.col(idx);
        p /= static_cast<double>(S.size());

        cov_out << "Facet " << f << " (" << S.size() << " pts): ";

        // Looping over scale factor
        for (int k = 1; k <= 10; ++k) 
        {
            double step = 0.1 * k;
            double x = std::pow(step, 1.0 / true_dim);

            // Local copy of polytope for each scaling
            HPoly P_loc = P;

            //Shifting and scaling
            P_loc.shift(p);
            const Eigen::MatrixXd T = (1.0 / x) * Eigen::MatrixXd::Identity(true_dim, true_dim);
            P_loc.linear_transformIt(T);

            //Parameters of new polytope
            const auto& A_sh = P_loc.get_mat();
            const auto& b_sh = P_loc.get_vec();

            // Points still in facet
            unsigned survivors = 0;
            for (int idx : S) {

                const Eigen::VectorXd q_shift = samples.col(idx) - p;
                bool inside = true;

                for (int j = 0; j < A_sh.rows(); ++j) {
                    if (j == f) continue;
                    if (A_sh.row(j).dot(q_shift) - b_sh[j] > tol) {
                        inside = false; break;
                    }
                }
                if (inside) ++survivors;
            }

            const double coverage = static_cast<double>(survivors) / S.size();
            cov_out << std::pow(x, true_dim) << ':' << coverage;
            if (k < 10) cov_out << ", ";
        }
        cov_out << '\n';  
    }
    cov_out.close();

    return 0;
}
