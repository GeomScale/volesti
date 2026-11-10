#include <ctime>
#include <iostream>
#include <list>
#include <vector>
#include <boost/random.hpp>
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/zpolytope.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sampling.hpp"
#include "generators/known_polytope_generators.h"
#include "generators/z_polytopes_generators.h"

typedef double NT;
typedef Cartesian<NT>::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

std::vector<unsigned> dims = {10, 20, 40, 60, 80, 100, 150, 200};

template <typename Polytope>
void run(const char *name, Polytope &P, unsigned d, unsigned calls = 500, unsigned samples = 20)
{
    boost::mt19937 gen(0);
    boost::normal_distribution<NT> N(0, 1);
    std::vector<Point> rs, vs;

    for (unsigned k = 0; k < calls; ++k) {
        Point r(d), v(d);
        for (unsigned i = 0; i < d; ++i) { 
            r.set_coord(i, N(gen));
            v.set_coord(i, N(gen));
        }
        rs.push_back(r*(0.25/(std::sqrt(NT(d))*r.length())));
        vs.push_back(v*(1.0/v.length()));
    }

    typename Polytope::VT dummy;

    std::clock_t t0 = std::clock();
    for (unsigned k = 0; k < calls; ++k) P.line_intersect(rs[k], vs[k]);
    std::clock_t t1 = std::clock();
    for (unsigned k = 0; k < calls; ++k) P.line_positive_intersect(rs[k], vs[k], dummy, dummy);
    std::clock_t t2 = std::clock();

    RNGType rng(d);
    rng.set_seed(0);
    std::list<Point> pts;
    uniform_sampling<BilliardWalk>(pts, P, rng, 1u, samples, Point(d), 0u);
    std::clock_t t3 = std::clock();

    auto us = [&](std::clock_t a, std::clock_t b) {return 1e6*(b-a)/CLOCKS_PER_SEC/calls; };

    std::cout << name << " d=" << d
              << " line_intersect=" << us(t0, t1) << "us"
              << " line_positive_intersect=" << us(t1, t2) << "us"
              << " billiard=" << double(t3-t2) / CLOCKS_PER_SEC << "s" 
              << std::endl;
}

int main()
{
    for (unsigned d : dims) {
        auto P = generate_cross<VPolytope<Point>>(d, true);
        run("cross", P, d);
    }

    for (unsigned d : dims) {
        auto Z = gen_zonotope_uniform<Zonotope<Point>, boost::mt19937>(d, 2*d, 0);
        run("zono", Z, d);
    }
    return 0;
}