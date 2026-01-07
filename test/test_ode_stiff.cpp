// Regression test for oscillatory ODE stability
// Reproduces issue#120 and validates fix in PR#375


#include "doctest.h"
#include <iostream>   // REQUIRED for runge_kutta.hpp
#include <cmath>
#include <vector>

#include "Eigen/Eigen"

#include "ode_solvers/ode_solvers.hpp"
#include "generators/known_polytope_generators.h"

TEST_CASE("rk4_oscillator_stability_regression") {

  using NT = double;
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef std::vector<Point> pts;
  typedef HPolytope<Point> Hpolytope;
  typedef std::vector<Hpolytope*> bounds;

  // x'' = -x  (harmonic oscillator)
  IsotropicQuadraticFunctor::parameters<NT> params;
  params.order = 2;
  params.alpha = 1;

  typedef IsotropicQuadraticFunctor::GradientFunctor<Point> func;
  func F(params);

  unsigned int dim = 1;

  // Initial conditions: x(0)=0, v(0)=1
  Point x0(dim);
  Point v0 = Point::all_ones(dim);

  pts state{x0, v0};

  RKODESolver<Point, NT, Hpolytope, func> solver(
      0.0,                 // initial time
      0.05,                // step size
      state,
      F,
      bounds{nullptr, nullptr}
  );

  const int steps = 2000;
  const NT amplitude_ub = 1.5;  // exact is 1, keep margin

  for (int i = 0; i < steps; ++i) {
    solver.step(i, true);

    NT x_norm = std::sqrt(solver.xs[0].dot(solver.xs[0]));
    CHECK(x_norm < amplitude_ub);
  }
}
