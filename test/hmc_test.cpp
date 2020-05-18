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
#include "euler.h"
#include "leapfrog.h"
#include "hmc.h"

template <typename NT>
void test_hmc_isotropic_gaussian(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef EulerODESolver<Point, NT, Vpolytope> Solver;
    // Isotropic gaussian
    NT L = 3;
    NT m = L;
    func neg_grad_f = [](pts x, NT t) { return (-1.0 / (2 * 3 * 3)) * x[0]; };
    std::function<NT(Point)> f = [](Point x) { return (0.5 / (2 * 3 * 3)) * x.dot(x); };
    unsigned int dim = 2;
    Point x0 = get_direction<RNGType, Point, NT>(dim, false);
    x0 = (1.0 / sqrt(3)) * x0;

    HMCSampler<Point, NT, RNGType, Vpolytope, Solver> hmc(neg_grad_f, f, x0, m, L, 0.1, 0.01, NULL);

    hmc.mix();

    for (int i = 0; i < 1000; i++) {
      Point p = hmc.sample();
      for (int j = 0; j < p.dimension(); j++) std::cout << p[j] << " ";
      std::cout << std::endl;
    }

}

template <typename NT>
void test_hmc_isotropic_gaussian_leapfrog(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef LeapfrogODESolver<Point, NT, Vpolytope> Solver;
    // Isotropic gaussian
    NT L = 3;
    NT m = L;
    func neg_grad_f = [](pts x, NT t) { return (-1.0 / (2 * 3 * 3)) * x[0]; };
    std::function<NT(Point)> f = [](Point x) { return (0.5 / (2 * 3 * 3)) * x.dot(x); };
    unsigned int dim = 2;
    Point x0 = get_direction<RNGType, Point, NT>(dim, false);
    x0 = (1.0 / sqrt(3)) * x0;

    HMCSampler<Point, NT, RNGType, Vpolytope, Solver> hmc(neg_grad_f, f, x0, m, L, 0.1, 0.01, NULL);

    for (int i = 0; i < 10000; i++) {
      Point p = hmc.sample();
      for (int j = 0; j < p.dimension(); j++) std::cout << p[j] << " ";
      std::cout << std::endl;
    }

}

template <typename NT>
void test_hmc_isotropic_truncated_gaussian(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef LeapfrogODESolver<Point, NT, Vpolytope> Solver;
    // Isotropic gaussian
    NT L = 1;
    NT m = L;
    func neg_grad_f = [](pts x, NT t) { return (-1.0) * x[0]; };
    std::function<NT(Point)> f = [](Point x) { return (0.5) * x.dot(x); };
    unsigned int dim = 2;
    Vpolytope P = gen_cube<Vpolytope>(dim, true);
    Point x0(dim);

    HMCSampler<Point, NT, RNGType, Vpolytope, Solver> hmc(neg_grad_f, f, x0, L, m, 0.01, 0.001, &P);

    hmc.samples(1000);

    for (int i = 0; i < 1000; i++) {
      Point p = hmc.sample();
      for (int j = 0; j < p.dimension(); j++) std::cout << p[j] << " ";
      std::cout << std::endl;
    }

}

template <typename NT>
void call_test_hmc() {
  test_hmc_isotropic_gaussian<double>();
  // test_hmc_isotropic_truncated_gaussian<double>();
  // test_hmc_isotropic_gaussian_leapfrog<double>();
}

TEST_CASE("hmc") {
  call_test_hmc<double>();
}
