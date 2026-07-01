// VolEsti (volume computation and sampling library)

// Copyright (c) 2025 Vissarion Fisikopoulos
// Copyright (c) 2025 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include <cmath>
#include <iostream>
#include <vector>
#include <functional>

#include "doctest.h"

#include "Eigen/Eigen"

#include "convex_bodies/convex_body.h"
#include "ode_solvers/generalized_leapfrog.hpp"
#include "ode_solvers/oracle_functors.hpp"
#include "cartesian_geom/cartesian_kernel.h"

// -------------------------------------------------------------------
// A harmonic-oscillator force functor that works for any number of
// position-velocity pairs in the coupled state vector.
// F(v_index, xs, t) returns -alpha * xs[v_index - 1], i.e. the
// acceleration for the pair whose velocity lives at index v_index.
// -------------------------------------------------------------------
template <typename Point>
struct HarmonicForce {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    NT alpha;

    HarmonicForce(NT alpha_) : alpha(alpha_) {}

    Point operator()(unsigned int const& i, pts const& xs, NT const& /*t*/) const {
        // i is the velocity index (odd).  The position for this pair is xs[i-1].
        return (-alpha) * xs[i - 1];
    }
};

// -------------------------------------------------------------------
// 1.  Single pair, no boundaries: verify the velocity-Verlet update
//     and approximate energy conservation.
// -------------------------------------------------------------------
template <typename NT>
void test_single_pair_no_boundaries() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;

    // d^2x/dt^2 = -x   (unit harmonic oscillator)
    HarmonicForce<Point> F(NT(1));

    Point x0(1);
    x0.set_coord(0, 1);            // position = 1
    Point v0(1);                    // velocity = 0
    std::vector<Point> q{x0, v0};

    GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
        NT(0),          // initial time
        NT(0.1),        // step size
        q,
        F,
        Bounds{nullptr, nullptr}   // no boundary constraints
    );

    // --- initial state ---
    CHECK(solver.dim == 1);
    CHECK(solver.t == NT(0));
    CHECK(solver.eta == NT(0.1));
    CHECK(solver.xs.size() == 2);
    CHECK(solver.num_steps == 0);
    CHECK(solver.num_reflections == 0);
    CHECK(std::abs(solver.xs[0][0] - NT(1)) < NT(1e-12));
    CHECK(std::abs(solver.xs[1][0] - NT(0)) < NT(1e-12));

    // --- one step ---
    // velocity-Verlet with eta = 0.1:
    //   v'   = 0      + 0.05 * (-1)           = -0.05
    //   x'   = 1      + 0.10 * (-0.05)        =  0.995
    //   v''  = -0.05  + 0.05 * (-0.995)       = -0.09975
    solver.step(0, true);

    CHECK(std::abs(solver.xs[0][0] - NT(0.995))   < NT(1e-12));
    CHECK(std::abs(solver.xs[1][0] - NT(-0.09975)) < NT(1e-12));
    CHECK(solver.t == NT(0.1));
    CHECK(solver.num_steps == 1);
    CHECK(solver.num_reflections == 0);

    // Energy H = v^2/2 + x^2/2  (harmonic oscillator Hamiltonian).
    // Velocity-Verlet is symplectic: H should be near-invariant.
    NT H_initial = NT(0.5);  // 0^2/2 + 1^2/2
    NT H_after   = solver.xs[1][0] * solver.xs[1][0] / NT(2)
                 + solver.xs[0][0] * solver.xs[0][0] / NT(2);
    CHECK(std::abs(H_after - H_initial) < NT(5e-4));

    // --- many steps: energy should not drift ---
    for (int i = 1; i < 100; ++i) {
        solver.step(i, true);
    }
    NT H_final = solver.xs[1][0] * solver.xs[1][0] / NT(2)
               + solver.xs[0][0] * solver.xs[0][0] / NT(2);
    CHECK(std::abs(H_final - H_initial) < NT(1e-2));
    CHECK(solver.num_steps == 100);
}

// -------------------------------------------------------------------
// 2.  Two independent pairs in the same solver.  Each pair is a 1-D
//     harmonic oscillator with a different initial condition.
// -------------------------------------------------------------------
template <typename NT>
void test_two_pairs() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;

    HarmonicForce<Point> F(NT(1));

    // pair 1: x = 1, v = 0
    // pair 2: x = 0.5, v = 0.5
    Point x1(1);  x1.set_coord(0, NT(1));
    Point v1(1);  // zero
    Point x2(1);  x2.set_coord(0, NT(0.5));
    Point v2(1);  v2.set_coord(0, NT(0.5));
    std::vector<Point> q{x1, v1, x2, v2};

    GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
        NT(0),
        NT(0.1),
        q,
        F,
        Bounds{nullptr, nullptr, nullptr, nullptr}
    );

    CHECK(solver.xs.size() == 4);
    CHECK(solver.dim == 1);

    solver.step(0, true);

    // pair 1: same as the single-pair test
    CHECK(std::abs(solver.xs[0][0] - NT(0.995))     < NT(1e-12));
    CHECK(std::abs(solver.xs[1][0] - NT(-0.09975))   < NT(1e-12));

    // pair 2: started at x=0.5, v=0.5
    //   v'   = 0.5     + 0.05 * (-0.5)  = 0.475
    //   x'   = 0.5     + 0.10 * 0.475   = 0.5475
    //   v''  = 0.475   + 0.05 * (-0.5475) = 0.447625
    CHECK(std::abs(solver.xs[2][0] - NT(0.5475))    < NT(1e-12));
    CHECK(std::abs(solver.xs[3][0] - NT(0.447625))   < NT(1e-12));
}

// -------------------------------------------------------------------
// 3.  Single pair constrained inside a unit ball (||x|| <= 1).
//     Position starts at the centre with a large velocity so that
//     the boundary is hit within the first step.
// -------------------------------------------------------------------
template <typename NT>
void test_constrained() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using NTp       = typename Point::FT;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;

    // Unit ball in 1-D:
    //   g(x) = x^2 - 1 <= 0  (interior)
    //   grad_g(x) = 2 * x
    std::vector<std::function<NTp(Point const&)>> gs;
    gs.push_back([](Point const& p) -> NTp {
        return p.dot(p) - NTp(1);
    });
    std::vector<std::function<Point(Point const&)>> grad_gs;
    grad_gs.push_back([](Point const& p) -> Point {
        return NTp(2) * p;
    });

    Convexbody ball(gs, grad_gs, 1);

    HarmonicForce<Point> F(NT(1));

    Point x0(1);                // position at centre (inside)
    Point v0(1);
    v0.set_coord(0, NT(20));    // large velocity toward the right boundary

    GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
        NT(0),
        NT(0.1),       // step size
        {x0, v0},
        F,
        Bounds{&ball, nullptr}
    );

    // Run enough steps to experience reflections and check that the
    // position never leaves the ball.
    unsigned int const STEPS = 100;
    for (unsigned int i = 0; i < STEPS; ++i) {
        solver.step(static_cast<int>(i), true);
        // The absolute position must stay inside the unit ball
        // (with a tiny tolerance for floating-point boundary rounding).
        CHECK(std::abs(solver.xs[0][0]) <= NTp(1) + NTp(1e-10));
    }

    CHECK(solver.num_reflections > 0);
}

// -------------------------------------------------------------------
// 4.  Adaptive stepping toggle and effect.
// -------------------------------------------------------------------
template <typename NT>
void test_adaptive_flag() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;

    HarmonicForce<Point> F(NT(1));

    Point x0(1);  x0.set_coord(0, NT(1));
    Point v0(1);
    std::vector<Point> q{x0, v0};

    {
        // Solver created with adaptive = true (the default)
        GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
            NT(0), NT(0.1), q, F, Bounds{nullptr, nullptr}
        );
        CHECK(solver.adaptive == true);
        solver.disable_adaptive();
        CHECK(solver.adaptive == false);
        solver.enable_adaptive();
        CHECK(solver.adaptive == true);
    }
    {
        // Solver created with adaptive = false
        GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
            NT(0), NT(0.1), q, F, Bounds{nullptr, nullptr}, false
        );
        CHECK(solver.adaptive == false);
    }
}

// -------------------------------------------------------------------
// 5.  State getters and setters.
// -------------------------------------------------------------------
template <typename NT>
void test_get_set_state() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;

    HarmonicForce<Point> F(NT(1));

    Point x0(1);  x0.set_coord(0, NT(1));
    Point v0(1);
    std::vector<Point> q{x0, v0};

    GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
        NT(0), NT(0.1), q, F, Bounds{nullptr, nullptr}
    );

    // get_state before any step
    Point p0 = solver.get_state(0);
    Point p1 = solver.get_state(1);
    CHECK(std::abs(p0[0] - NT(1)) < NT(1e-12));
    CHECK(std::abs(p1[0] - NT(0)) < NT(1e-12));

    // set_state
    Point new_pos(1);  new_pos.set_coord(0, NT(0.5));
    Point new_vel(1);  new_vel.set_coord(0, NT(1.5));
    solver.set_state(0, new_pos);
    solver.set_state(1, new_vel);
    CHECK(std::abs(solver.xs[0][0] - NT(0.5)) < NT(1e-12));
    CHECK(std::abs(solver.xs[1][0] - NT(1.5)) < NT(1e-12));
}

// -------------------------------------------------------------------
// 6.  Using the library-provided IsotropicQuadraticFunctor with
//     order=2 (second-order ODE, i.e. d^2x/dt^2 = -alpha * x).
// -------------------------------------------------------------------
template <typename NT>
void test_isotropic_quadratic_functor() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;
    using GradFunc  = IsotropicQuadraticFunctor::GradientFunctor<Point>;
    using FuncParams = IsotropicQuadraticFunctor::parameters<NT>;

    FuncParams params(NT(1.0), 2);  // alpha = 1, order = 2
    GradFunc F(params);

    Point x0(1);  x0.set_coord(0, NT(1));
    Point v0(1);
    std::vector<Point> q{x0, v0};

    GeneralizedLeapfrogODESolver<Point, NT, Convexbody, GradFunc> solver(
        NT(0), NT(0.1), q, F, Bounds{nullptr, nullptr}
    );

    CHECK(solver.dim == 1);
    CHECK(solver.eta == NT(0.1));

    // single step (same expected values as test_single_pair_no_boundaries)
    solver.step(0, true);
    CHECK(std::abs(solver.xs[0][0] - NT(0.995))   < NT(1e-12));
    CHECK(std::abs(solver.xs[1][0] - NT(-0.09975)) < NT(1e-12));
}

// -------------------------------------------------------------------
// 7.  Two-dimensional harmonic oscillator without boundaries.
// -------------------------------------------------------------------
template <typename NT>
void test_2d_no_boundaries() {
    using Kernel    = Cartesian<NT>;
    using Point     = typename Kernel::Point;
    using Convexbody = ConvexBody<Point>;
    using Bounds    = std::vector<Convexbody*>;

    HarmonicForce<Point> F(NT(1));

    Point x0(2);
    x0.set_coord(0, NT(1));
    x0.set_coord(1, NT(-1));

    Point v0(2);
    v0.set_coord(0, NT(0.5));
    v0.set_coord(1, NT(0.5));

    std::vector<Point> q{x0, v0};

    GeneralizedLeapfrogODESolver<Point, NT, Convexbody, HarmonicForce<Point>> solver(
        NT(0), NT(0.05), q, F, Bounds{nullptr, nullptr}
    );

    CHECK(solver.dim == 2);
    CHECK(solver.xs.size() == 2);

    // The functor returns acceleration = -position.  With a small step
    // the norm of the coupled state should stay near its initial value.
    NT H0 = x0.dot(x0) / NT(2) + v0.dot(v0) / NT(2);  // 1^2+(-1)^2 = 2, 0.5^2+0.5^2=0.5 => H0 = 1 + 0.25 = 1.25

    for (int i = 0; i < 50; ++i) {
        solver.step(i, true);
    }

    NT H = solver.xs[0].dot(solver.xs[0]) / NT(2)
         + solver.xs[1].dot(solver.xs[1]) / NT(2);
    CHECK(std::abs(H - H0) / H0 < NT(1e-2));
}

// ===================================================================
//  TEST CASE REGISTRATIONS
// ===================================================================

TEST_CASE("generalized_leapfrog_no_boundaries") {
    std::cout << "--- Single pair, no boundaries (1-D) ---" << std::endl;
    test_single_pair_no_boundaries<double>();
}

TEST_CASE("generalized_leapfrog_two_pairs") {
    std::cout << "--- Two independent pairs ---" << std::endl;
    test_two_pairs<double>();
}

TEST_CASE("generalized_leapfrog_constrained") {
    std::cout << "--- Constrained (unit ball) ---" << std::endl;
    test_constrained<double>();
}

TEST_CASE("generalized_leapfrog_adaptive_flag") {
    std::cout << "--- Adaptive flag toggle ---" << std::endl;
    test_adaptive_flag<double>();
}

TEST_CASE("generalized_leapfrog_get_set_state") {
    std::cout << "--- State getters / setters ---" << std::endl;
    test_get_set_state<double>();
}

TEST_CASE("generalized_leapfrog_isotropic_quadratic_functor") {
    std::cout << "--- IsotropicQuadraticFunctor ---" << std::endl;
    test_isotropic_quadratic_functor<double>();
}

TEST_CASE("generalized_leapfrog_2d") {
    std::cout << "--- 2-D harmonic oscillator ---" << std::endl;
    test_2d_no_boundaries<double>();
}
