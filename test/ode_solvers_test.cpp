#include "ode_solvers.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include "Eigen/Eigen"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "known_polytope_generators.h"
#include <string>
#include <typeinfo>
#include "samplers.h"
#include "doctest.h"

template <typename NT>
void test_euler(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    funcs Fs;

    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);

    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Vpolytope> euler_solver = EulerODESolver<Point, NT, Vpolytope>(0, 0.01, q, Fs);
    euler_solver.steps(1000);
    NT err=0.001;
    NT error = std::abs(euler_solver.xs[0][0]);
    CHECK(error < err);
}

template <typename NT>
void call_test_euler() {
  std::cout << "--- Testing solution to dx / dt = -x" << std::endl;
  test_euler<NT>();
}

TEST_CASE("euler") {
  call_test_euler<double>();
}
