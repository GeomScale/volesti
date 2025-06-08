<<<<<<< Updated upstream
/*  shakeandbake_examples_modes_square.cpp
 *  Shake-and-Bake za kvadrat definisan na slici:
 *    – vrhovi kvadrata su: (-1,-1), (1,-1), (1,1), (-1,1)
 *    – ograničenja (H-poliedar):
 *        x ≤  1        (desna strana)
 *       −x ≤  1        (leva strana)
 *        y ≤  1        (gornja strana)
 *       −y ≤  1        (donja strana)
 *  – Svaki mod startuje NA FASETI, ali volesti preskače λ=0 za tu fasetu.
 *  – Izlazi: sb_square_orig.txt, sb_square_limp.txt, sb_square_run.txt.
 */

#include <iostream>
#include <fstream>
#include <chrono>
=======
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
>>>>>>> Stashed changes
#include <Eigen/Eigen>
#include <boost/random.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "random_walks/shake_and_bake_walk.hpp"
<<<<<<< Updated upstream
=======
#include "generators/known_polytope_generators.h" 
>>>>>>> Stashed changes

using NT     = double;
using Kernel = Cartesian<NT>;
using Point  = Kernel::Point;
<<<<<<< Updated upstream
using MT     = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
using VT     = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
using RNG    = BoostRandomNumberGenerator<boost::random::mt19937, NT>;
using HPoly  = HPolytope<Point>;

/**
 *  Pravi H-poliedar koji opisuje kvadrat sa temenima:
 *    A = (-1,-1), B = (1,-1), C = (1,1), D = (-1,1)
 *
 *  Ograničenja su:
 *    (1)  x ≤  1      → [1  0]·[x y]^T ≤ 1
 *    (2) −x ≤  1      → [−1  0]·[x y]^T ≤ 1
 *    (3)  y ≤  1      → [0  1]·[x y]^T ≤ 1
 *    (4) −y ≤  1      → [0 −1]·[x y]^T ≤ 1
 */
static HPoly square_polytope()
{
    MT A(4, 2);
    A <<   1,  0,   // x ≤  1
          -1,  0,   // −x ≤  1
           0,  1,   // y ≤  1
           0, -1;   // −y ≤  1

    VT b(4);
    b << 1, 1, 1, 1;

    return HPoly(2, A, b);
}

int main()
{
    HPoly P = square_polytope();
    RNG rng(P.dimension());

    // Početna tačka na gornjoj strani kvadrata (y = 1) – uzimamo sredinu: (0,1)
    Point p0(P.dimension());
    p0.set_coord(0, 0.0);  // x = 0
    p0.set_coord(1, 1.0);  // y = 1

    const unsigned walk_len  = 10;
    const unsigned n_samples = 100;

    struct ModeInfo {
        SBWalk::Mode mode;
        const char*  filename;
    } modes[] = {
        { SBWalk::Mode::Original, "sb_square_orig.txt" },
        { SBWalk::Mode::Limping,  "sb_square_limp.txt"  },
        { SBWalk::Mode::Running,  "sb_square_run.txt"   }
    };

    for (auto &mi : modes) {
        std::ofstream outfile(mi.filename);
        if (!outfile.is_open()) {
            std::cerr << "Greška: ne mogu da otvorim \"" << mi.filename << "\"\n";
            continue;
        }

        // Inicijalizujemo shake-and-bake walker u datom modu
        using Walker = SBWalk::Walk<HPoly, RNG>;
        Walker walk(P, p0, rng, mi.mode);

        // Zabeležimo trenutak pre generisanja uzoraka
        auto t0 = std::chrono::high_resolution_clock::now();

        Point d1(P.dimension());
        for (unsigned i = 0; i < n_samples; ++i) {
            // Izvrši walk_len koraka i uzmi uzorak
            walk.apply(P, d1, walk_len, rng);
            const Point& q = walk.getCurrentPoint();
            outfile << q[0] << ' ' << q[1] << "\n";
        }

        // Zabeležimo trenutak posle generisanja
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t1 - t0;

        outfile.close();

        // Ispišemo na konzolu podatke o modu i koliko je trajalo (u sekundama)
        std::cout << "Mod ";
        switch (mi.mode) {
            case SBWalk::Mode::Original: std::cout << "Original"; break;
            case SBWalk::Mode::Limping:  std::cout << "Limping";  break;
            case SBWalk::Mode::Running:  std::cout << "Running";  break;
            default:                     std::cout << "Nepoznat"; break;
        }
        std::cout << ": " << n_samples << " uzoraka u \"" << mi.filename 
                  << "\" generisano za " << elapsed.count() << " s\n";
=======
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

    // Generating polytope 
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

    // Calculating walk length and number of samples according to number of dimensions
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
        { SBWalk::Mode::Limping,  "_limp.txt" },
        { SBWalk::Mode::Running,  "_run.txt"  }
    };

    for (auto &mi : modes) {
        std::string fname = base + mi.suffix;
        std::ofstream outfile(fname);
        if (!outfile.is_open()) {
            std::cerr << "Greška: ne mogu da otvorim \"" << fname << "\"\n";
            continue;
        }

        using Walker = SBWalk::Walk<HPoly, RNG>;
        Walker walk(P, p0, rng, mi.mode);

        auto t0 = std::chrono::high_resolution_clock::now();
        Point d1(dim);

        for (unsigned i = 0; i < n_samples; ++i) {
            walk.apply(P, d1, walk_len, rng);
            const Point& q = walk.getCurrentPoint();
            for (unsigned j = 0; j < dim; ++j)
                outfile << q[j] << (j + 1 < dim ? ' ' : '\n');
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        outfile.close();

        std::cout << "Mode ";
        switch (mi.mode) {
            case SBWalk::Mode::Original: std::cout << " Shake and Bake Original"; break;
            case SBWalk::Mode::Limping:  std::cout << " Shake and Bake Limping";  break;
            case SBWalk::Mode::Running:  std::cout << " Shake and Bake Running";  break;
            default:                     std::cout << "Unknown"; break;
        }
        std::cout << ": " << n_samples << " samples (walk_len=" << walk_len 
                  << ") u \"" << fname << "\" generated for " << elapsed << " s\n";
>>>>>>> Stashed changes
    }

    return 0;
}
