// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

// The interrupt hook VOLESTI_CHECK_INTERRUPT, defined here as a binding
// would define it: before any volesti header. It counts the calls and
// "interrupts" (throws) after a limit.

#include <list>

struct interrupted {};
static unsigned long hook_calls = 0;
static unsigned long hook_limit = 0;
#define VOLESTI_CHECK_INTERRUPT() \
    do { if (++hook_calls > hook_limit) throw interrupted(); } while (0)

#include "doctest.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "generators/known_polytope_generators.h"
#include "sampling/sampling.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

template <typename WalkType>
void sample_cube(unsigned int walk_len, unsigned int rnum)
{
    Hpolytope P = generate_cube<Hpolytope>(10, false);
    P.ComputeInnerBall();
    RNGType rng(P.dimension());
    std::list<Point> randPoints;
    uniform_sampling<WalkType>(randPoints, P, rng, walk_len, rnum,
                               Point(P.dimension()), 0);
}

template <typename WalkType>
void call_test_interrupt()
{
    // A walk of 10^9 steps per point stops at the first throwing check
    hook_calls = 0;
    hook_limit = 1000;
    CHECK_THROWS_AS(sample_cube<WalkType>(1000000000, 1), interrupted);
    CHECK(hook_calls == hook_limit + 1);
}

TEST_CASE("hook_called_once_per_step") {
    hook_calls = 0;
    hook_limit = std::numeric_limits<unsigned long>::max();
    sample_cube<CDHRWalk>(5, 100);
    CHECK(hook_calls == 5 * 100);
}

TEST_CASE("interrupt") {
    call_test_interrupt<CDHRWalk>();
    call_test_interrupt<RDHRWalk>();
    call_test_interrupt<BallWalk>();
    call_test_interrupt<BilliardWalk>();
    call_test_interrupt<AcceleratedBilliardWalk>();
}
