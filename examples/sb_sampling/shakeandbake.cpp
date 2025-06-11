
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <Eigen/Eigen>
#include <boost/random.hpp>

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

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Template: " << argv[0] << " <cube|simplex> <number of dimensions>\n";
        return 1;
    }
    std::string shape = argv[1];
    unsigned    dim   = std::stoi(argv[2]);

    HPoly P;
    if (shape == "cube") {
        P = generate_cube<HPoly>(dim, false);
    }
    else if (shape == "simplex") {
        P = generate_simplex<HPoly>(dim, false);
    }
    else {
        std::cerr << "Unknown polytope '" << shape << "'. Use 'cube' or 'simplex'.\n";
        return 1;
    }

    unsigned walk_len  = 2 * dim * dim;
    unsigned n_samples = 100 * dim;

    RNG rng(P.dimension());
    Point p0(dim);
    if (shape == "cube") {
        for (unsigned i = 0; i < dim-1; ++i) p0.set_coord(i, 0.0);
        p0.set_coord(dim-1, 1.0);
    } else {
        for (unsigned i = 0; i < dim; ++i)
            p0.set_coord(i, -1.0/static_cast<NT>(dim));
        p0.set_coord(0, 0.0);
    }

    std::string base = "sb_" + shape + "_" + std::to_string(dim);
    struct ModeInfo {
        SBWalk::Mode mode;
        const char*  suffix;
    } modes[] = {
        { SBWalk::Mode::Original, "_orig.txt" },
        { SBWalk::Mode::Limping,  "_limp.txt"  },
        { SBWalk::Mode::Running,  "_run.txt"   }
    };

    for (auto &mi : modes) {
        std::string fname = base + mi.suffix;
        std::ofstream outfile(fname);

        using Walker = SBWalk::Walk<HPoly, RNG>;
        Walker walk(P, p0, rng, mi.mode);

        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> samples(dim, n_samples);

        auto t0 = std::chrono::high_resolution_clock::now();
        Point d1(dim);

        for (unsigned i = 0; i < n_samples; ++i) {
            walk.apply(P, d1, walk_len, rng);
            const Point& q = walk.getCurrentPoint();
            for (unsigned j = 0; j < dim; ++j) {
                outfile << q[j] << (j + 1 < dim ? ' ' : '\n');
                samples(j, i) = q[j];
            }
        }
        outfile.close();

        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();

        std::cout << "Mode ";
        switch (mi.mode) {
            case SBWalk::Mode::Original: std::cout << "Shake and Bake Original"; break;
            case SBWalk::Mode::Limping:  std::cout << "Shake and Bake Limping";  break;
            case SBWalk::Mode::Running:  std::cout << "Shake and Bake Running";  break;
            default:                     std::cout << "Unknown";             break;
        }
        std::cout << ": " << n_samples << " samples (walk_len=" << walk_len 
                  << ") generated in " << elapsed << " s\n";

        unsigned int min_ess = 0;
        std::cout << "== Diagnostics per dimension ==\n";
        print_diagnostics<NT, Eigen::VectorXd, decltype(samples), std::ostream>(samples, min_ess, std::cout);
        std::cout << "Minimum ESS (for thinning): " << min_ess << "\n\n";
    }

    return 0;
}
