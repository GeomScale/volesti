#include <iostream>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "ellipsoid.h"
#include "ellipsoidintersectconvex.h"
#include "sampling/ellipsoid.hpp"
#include "random_walks/random_walks.hpp"

template <typename NT>
int run_ellipsoid_sampling() {
    typedef Cartesian<NT>    Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point    Point;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    unsigned int dim = 2;
    MT A(2, 2);
    A << 0.25, 0.75,
         0.75, 3.25;

    Ellipsoid<Point> ell(A);    // origin centered ellipsoid
    int num_points = 5000;
    Point p(dim);
    RNGType rng(dim);

    std::cout << "Sampling from ellipsoid..." << std::endl;

    for (int i=0; i<num_points; ++i) {
        p = GetPointInDellipsoid<Point>::apply(dim, ell, rng);
        p.print();
    }

    std::cout << "Sampling from ellipsoid completed." << std::endl;

    return 0;
}

template <typename NT>
int run_ellipsoid_polytope_sampling() {
    typedef Cartesian<NT>    Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point    Point;
    typedef HPolytope <Point> Hpolytope;
    typedef Ellipsoid <Point> CEllipsoid;
    typedef EllipsoidIntersectPolytope <Hpolytope, CEllipsoid> EllipsoidIntPolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef AcceleratedBilliardWalk::template Walk<EllipsoidIntPolytope, RNGType> AcceleratedBilliardWalkType;

    PushBackWalkPolicy push_back_policy;
    unsigned int dim = 2;
    MT E(2, 2);
    E << 0.25, 0.75,
         0.75, 3.25;

    Ellipsoid<Point> ell(E);    // origin centered ellipsoid

    MT A(4,2);
    A <<  1.0,  0.0,
          0.0,  1.0,
         -1.0,  0.0,
          0.0, -1.0;
    VT b = VT::Ones(4);
    Hpolytope P(dim, A, b);

    EllipsoidIntPolytope EP(P, ell);

    int num_points = 5000;
    RNGType rng(dim);
    Point p(dim);

    std::cout << "Sampling from the intersection between ellipsoid and polytope..." << std::endl;

    typedef RandomPointGenerator<AcceleratedBilliardWalkType> Generator;
    std::vector<Point> randPoints;
    Point q = Point(VT::Zero(dim)); // origin
    Generator::apply(EP, q, num_points, 1,
                     randPoints, push_back_policy, rng);

    for (int i=0; i<randPoints.size(); ++i) {
        p = randPoints[i];
        p.print();
    }

    std::cout << "Sampling from the intersection completed." << std::endl;

    return 0;
}

struct CustomFunctor {

  // Custom density with neg log prob equal to || x ||^2 + 1^T x
  template <
      typename NT
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number

    parameters() : order(2), L(2), m(2), kappa(1) {};

    parameters(unsigned int order_) :
      order(order),
      L(2),
      m(2),
      kappa(1)
    {}
  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> params;

    GradientFunctor() {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * Point::all_ones(xs[0].dimension());
        y = y + (-2.0) * xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> params;

    FunctionFunctor() {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      return x.dot(x) + x.sum();
    }

  };

};


template <typename NT>
int run_ellipsoid_polytope_nuts_sampling() {
    typedef Cartesian<NT>    Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef HPolytope <Point> Hpolytope;
    typedef Ellipsoid <Point> CEllipsoid;
    typedef EllipsoidIntersectPolytope <Hpolytope, CEllipsoid> EllipsoidIntPolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, EllipsoidIntPolytope, NegativeGradientFunctor> Solver;

    
    PushBackWalkPolicy push_back_policy;
    unsigned int dim = 2;
    MT E(2, 2);
    E << 0.25, 0.75,
         0.75, 3.25;

    Ellipsoid<Point> ell(E);    // origin centered ellipsoid

    MT A(4,2);
    A <<  1.0,  0.0,
          0.0,  1.0,
         -1.0,  0.0,
          0.0, -1.0;
    VT b = VT::Ones(4);

    Hpolytope P(dim, A, b);
    EllipsoidIntPolytope EP(P, ell);

    NegativeGradientFunctor F;
    NegativeLogprobFunctor f;
    NutsHamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
    
    std::cout << "Starting NUTS example" << std::endl;
    std::cout << "eta0 before burnin: " << hmc_params.eta << std::endl;
    
    int num_points = 5000;
    bool automatic_burnin = false;
    RNGType rng(dim);
    Point p(dim), q(dim);

    NutsHamiltonianMonteCarloWalk::Walk<Point, EllipsoidIntPolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
        hmc(&EP, p, F, f, hmc_params, automatic_burnin);
    hmc.burnin(rng);

    std::cout << "eta after burnin: " << hmc.get_eta_solver() << std::endl;
    std::cout << "Sampling with NUTS from the intersection between ellipsoid and polytope..." << std::endl;

    for (int i = 0; i < num_points; i++) {
        hmc.apply(rng);
        q = hmc.x.getCoefficients();
        q.print();
    }
    
    std::cout << "Sampling with NUTS from the intersection completed." << std::endl;

    return 0;
}


int main() {
  run_ellipsoid_sampling<double>();
  run_ellipsoid_polytope_sampling<double>();
  run_ellipsoid_polytope_nuts_sampling<double>();
  return 0;
}
