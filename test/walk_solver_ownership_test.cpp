// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

// The HMC, NUTS-HMC and Langevin walks each create an ODE/SDE solver when they
// are constructed. These tests construct a walk, use it and let it go out of
// scope. Run under LeakSanitizer (-fsanitize=address) they detect a solver that
// is not released together with its walk; in a regular build they only check that
// the walks run and produce finite points.

#include "doctest.h"

#include <cmath>
#include <vector>

#include <boost/random.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "generators/known_polytope_generators.h"
#include "ode_solvers/ode_solvers.hpp"
#include "random_walks/random_walks.hpp"
#include "random_walks/hamiltonian_monte_carlo_walk.hpp"
#include "random_walks/nuts_hmc_walk.hpp"
#include "random_walks/langevin_walk.hpp"

typedef Cartesian<double>                                          Kernel;
typedef Kernel::Point                                              Point;
typedef HPolytope<Point>                                           Hpolytope;
typedef BoostRandomNumberGenerator<boost::mt19937, double, 3>      RNG;

// Standard Gaussian density: negative log-probability f(x) = |x|^2 / 2
struct StandardGaussian {
    struct parameters {
        unsigned int order = 2;
        double L = 1.0;      // Lipschitz constant of the gradient
        double m = 1.0;      // strong convexity constant
        double kappa = 1.0;  // condition number
    };

    struct GradientFunctor {
        typedef std::vector<Point> pts;
        parameters params;
        // state vector index i: the derivative of x is v, the derivative of v is -grad f(x)
        Point operator()(unsigned int const& i, pts const& xs, double const&) const {
            if (i == params.order - 1) {
                return (-1.0) * xs[0];
            }
            return xs[i + 1];
        }
    };

    struct FunctionFunctor {
        parameters params;
        double operator()(Point const& x) const { return 0.5 * x.dot(x); }
    };
};

template <class Walk>
static void run_walk(Walk& walk, RNG& rng)
{
    for (int i = 0; i < 20; ++i) {
        walk.apply(rng, 1);
        for (unsigned int j = 0; j < walk.x.dimension(); ++j) {
            CHECK(std::isfinite(walk.x[j]));
        }
    }
}

TEST_CASE("hmc_walk_releases_solver") {
    typedef StandardGaussian::GradientFunctor Grad;
    typedef StandardGaussian::FunctionFunctor Func;
    typedef LeapfrogODESolver<Point, double, Hpolytope, Grad> Solver;

    unsigned int const dim = 3;
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Grad F; Func f; RNG rng(1); Point x0(dim);
    HamiltonianMonteCarloWalk::parameters<double, Grad> params(F, dim);
    HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RNG, Grad, Func, Solver>
        walk(&P, x0, F, f, params);
    run_walk(walk, rng);
}

TEST_CASE("nuts_walk_releases_solver") {
    typedef StandardGaussian::GradientFunctor Grad;
    typedef StandardGaussian::FunctionFunctor Func;
    typedef LeapfrogODESolver<Point, double, Hpolytope, Grad> Solver;

    unsigned int const dim = 3;
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Grad F; Func f; RNG rng(1); Point x0(dim);
    NutsHamiltonianMonteCarloWalk::parameters<double, Grad> params(F, dim);
    NutsHamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RNG, Grad, Func, Solver>
        walk(&P, x0, F, f, params, false);
    run_walk(walk, rng);
}

TEST_CASE("langevin_walk_releases_solver") {
    typedef StandardGaussian::GradientFunctor Grad;
    typedef StandardGaussian::FunctionFunctor Func;
    typedef RandomizedMipointSDESolver<Point, double, Hpolytope, Grad, RNG> Solver;

    unsigned int const dim = 3;
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Grad F; Func f; RNG rng(1); Point x0(dim);
    UnderdampedLangevinWalk::parameters<double, Grad> params(F, dim);
    UnderdampedLangevinWalk::Walk<Point, Hpolytope, RNG, Grad, Func, Solver>
        walk(&P, x0, F, f, params);
    run_walk(walk, rng);
}
