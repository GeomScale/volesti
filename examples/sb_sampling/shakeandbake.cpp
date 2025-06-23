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

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage : " << argv[0]
                  << " <cube|simplex|birkhoff|iSDY_1059> <dimension> [epsilon]\n";
        return 1;
    }

    std::string shape = argv[1];
    unsigned    cli_n = std::stoi(argv[2]);
    double eps = (argc > 3 ? std::stod(argv[3]) : 1e-7);

    HPoly P;
    if (shape == "cube")               P = generate_cube<HPoly>(cli_n, false);
    else if (shape == "simplex")       P = generate_simplex<HPoly>(cli_n, false);
    else if (shape == "birkhoff")      P = generate_birkhoff<HPoly>(cli_n);
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
    if (shape == "cube" || shape == "simplex") {
        walk_len      = 50;
        n_samples     = 100000;
        burn_in_iters = 0;
    } else if (shape == "birkhoff") {
       walk_len      = 20 * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 5  * true_dim;
    } else {
        walk_len      = 20 * true_dim;
        n_samples     = 100 * true_dim;
        burn_in_iters = 20 * true_dim;
    }

    std::cout << "Parameters: walk_len=" << walk_len
              << ", n_samples=" << n_samples
              << ", burn_in_iters=" << burn_in_iters
              << " (dim=" << true_dim << ") eps=" << eps << "\n";

    RNG rng(true_dim);
    std::string base = "sb_" + shape + "_" + std::to_string(cli_n);
    std::ofstream out(base + "_run.txt");

    using Walker = ShakeAndBakeWalk::Walk<HPoly, RNG>;
    Walker walk(P, rng, ShakeAndBakeWalk::Running);

    auto A_full = P.get_mat();
    auto b_full = P.get_vec();
    size_t m = A_full.rows();

    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples(true_dim, n_samples);
    std::vector<int> facet_id(n_samples, -1);

    // burn-in
    Point tmp0(true_dim);
    for (unsigned i = 0; i < burn_in_iters; ++i)
        walk.apply(P, tmp0, walk_len, rng);

    // sampling
    for (unsigned i = 0; i < n_samples; ++i) {
        Point tmp(true_dim);
        walk.apply(P, tmp, walk_len, rng);
        const Point &q = walk.getCurrentPoint();
        Eigen::VectorXd qv(true_dim);
        for (unsigned d=0; d<true_dim; ++d) qv[d] = q[d];
        auto Aq = A_full * qv;
        for (size_t k=0; k<m; ++k) {
            if (std::abs(Aq[k] - b_full[k]) < eps) {
                facet_id[i] = k;
                break;
            }
        }
        for (unsigned d=0; d<true_dim; ++d) {
            out << q[d] << (d+1<true_dim ? ' ' : '\n');
            samples(d,i) = q[d];
        }
    }
    out.close();

    std::cout << "Generated " << n_samples << " samples in "
              << walk_len << " steps each.\n";


    const double min_ratio = 0.01;
    std::ofstream cov_out("coverage.txt");

    cov_out    << "\nScaling coverage by facet (skip <" << min_ratio << "):\n";

    for (int f = 0; f < int(m); ++f)
    {
        //Samples S on facet f
        std::vector<int> S;
        for (unsigned i = 0; i < n_samples; ++i)
            if (facet_id[i] == f) S.push_back(i);

        //Setting threshold
        const double ratio = double(S.size()) / n_samples;
        if (ratio < min_ratio) {
            cov_out << "facet " << f << " skipped (" << ratio << ")\n";
            continue;
        }

        //Finding the mean of the samples
        Eigen::VectorXd p = Eigen::VectorXd::Zero(true_dim);
        for (int idx : S) p += samples.col(idx);
        p /= double(S.size());

        cov_out  << "Facet " << f << " (" << S.size() << " pts): ";

        const Eigen::VectorXd Af_orig = A_full.row(f);
        const double          bf_orig = b_full[f];

        // Loop over scale factors
        for (int step = 0; step <= 50; ++step)
        {
            const double x = 0.01 + 0.99 * step / 50.0;   // shrink factor 

            HPoly P_loc = P;          // copy for each scale 
            P_loc.shift(p);

            //Shrinking
            const Eigen::MatrixXd T = (1.0 / x) *
                                    Eigen::MatrixXd::Identity(true_dim, true_dim);
            P_loc.linear_transformIt(T);

            const Eigen::MatrixXd& A_sh = P_loc.get_mat();
            const Eigen::VectorXd& b_sh = P_loc.get_vec();
            const Eigen::VectorXd  Af   = A_sh.row(f);
            const double           bf   = b_sh[f];

            unsigned survivors = 0;

            for (int idx : S)
            {
                const Eigen::VectorXd q = samples.col(idx) - p;

                bool inside = true;
                for (int j = 0; j < A_sh.rows(); ++j)
                    if (A_sh.row(j).dot(q) - b_sh[j] > eps) { inside = false; break; }
                if (!inside) continue;

                if (std::abs(Af.dot(q) - bf) < eps) ++survivors;
            }

            const double coverage = double(survivors) / S.size();
            cov_out  << x << ':' << coverage << (step < 50 ? ", " : "\n");
        }
    }
    cov_out.close();

    return 0;
}
