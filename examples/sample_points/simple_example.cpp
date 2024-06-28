#include <chrono>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"

#include "diagnostics/diagnostics.hpp"
//todo make headers automomous
//#include "diagnostics/effective_sample_size.hpp"
//#include "diagnostics/print_diagnostics.hpp"

#include "generators/known_polytope_generators.h"

#include "random_walks/random_walks.hpp"

#include "sampling/sample_points.hpp"
#include "sampling/sampling.hpp"

#include "volume/volume_cooling_balls.hpp"

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
typedef HPolytope<Point> HPolytopeType;

using NT = double;
using MT = Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>;
using VT = Eigen::Matrix<NT,Eigen::Dynamic,1>;

template <typename PolytopeOrProblem, typename Walk, typename Distribution>
void sample_points_eigen_matrix(PolytopeOrProblem const& HP, Point const& q, Walk const& walk,
                                Distribution const& distr, RNGType rng, int walk_len, int rnum,
                                int nburns)
{
    MT samples(q.dimension(), rnum);

    auto start = std::chrono::steady_clock::now();

    sample_points(HP, q, walk, distr, rng, walk_len, rnum, nburns, samples);

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << " time= " << time;

    // sample stats
    unsigned int min_ess;
    auto score = effective_sample_size<NT, VT, MT>(samples, min_ess);
    std::cout << " ess= "  << min_ess << std::endl;
    //print_diagnostics<NT, VT, MT>(samples, min_ess, std::cerr);
}

struct CustomFunctor {

  // Custom density with neg log prob equal to c^T x
  template
  <
      typename NT,
      typename Point
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number
    Point x0;

    parameters(Point x0_) : order(2), L(1), m(1), kappa(1), x0(x0_) {};

  };

  template <typename Point>
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * (xs[0] - params.x0);
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }
  };

  template<typename Point>
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      Point y = x - params.x0;
      return 0.5 * y.dot(y);
    }
  };
};

inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

template< typename SpMat, typename VT>
void load_crhmc_problem(SpMat &A, VT &b, VT &lb, VT &ub, int &dimension,
                        std::string problem_name) {
   {
    std::string fileName("../crhmc_sampling/data/");
    fileName.append(problem_name);
    fileName.append(".mm");
    if(!exists_check(fileName)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat X;
    loadMarket(X, fileName);
    int m = X.rows();
    dimension = X.cols() - 1;
    A = X.leftCols(dimension);
    b = VT(X.col(dimension));
  }
  {
    std::string fileName("../crhmc_sampling/data/");
    fileName.append(problem_name);
    fileName.append("_bounds.mm");
    if(!exists_check(fileName)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat bounds;
    loadMarket(bounds, fileName);
    lb = VT(bounds.col(0));
    ub = VT(bounds.col(1));
  }
}


int main() {
    // NEW INTERFACE Sampling

    // Inputs:

    // Generating a 3-dimensional cube centered at origin
    HPolytopeType HP = generate_cube<HPolytopeType>(10, false);
    std::cout<<"Polytope: \n";
    HP.ComputeInnerBall();
    //HP.print();
    //std::cout<<"\n";

    // Setup parameters for sampling
    Point q(HP.dimension());
    RNGType rng(HP.dimension());

    // Generating a sparse polytope/problem
    using SpMat = Eigen::SparseMatrix<NT>;
    using ConstraintProblem =constraint_problem<SpMat, Point>;
    std::string problem_name("simplex3");
    std::cerr << "CRHMC on " << problem_name << "\n";
    SpMat As;
    VT b, lb, ub;
    int dimension;
    load_crhmc_problem(As, b, lb, ub, dimension, problem_name);
    ConstraintProblem problem = ConstraintProblem(dimension);
    problem.set_equality_constraints(As, b);
    problem.set_bounds(lb, ub);


    // Walks
    AcceleratedBilliardWalk abill_walk;
    AcceleratedBilliardWalk abill_walk_custom(10); //user defined walk parameters
    BallWalk ball_walk;
    BilliardWalk bill_walk;
    CDHRWalk cdhr_walk;
    DikinWalk dikin_walk;
    JohnWalk john_walk;
    RDHRWalk rdhr_walk;
    VaidyaWalk vaidya_walk;

    GaussianBallWalk gball_walk;
    GaussianCDHRWalk gcdhr_walk;
    GaussianRDHRWalk grdhr_walk;
    GaussianHamiltonianMonteCarloExactWalk ghmc_walk;

    GaussianAcceleratedBilliardWalk gbill_walk;

    ExponentialHamiltonianMonteCarloExactWalk ehmc_walk;

    HamiltonianMonteCarloWalk hmc_walk;
    NutsHamiltonianMonteCarloWalk nhmc_walk;
    UnderdampedLangevinWalk uld_walk;
    CRHMCWalk crhmc_walk;

    // Distributions

    // 1. Uniform
    UniformDistribution udistr{};

    // 2. Spherical
    SphericalGaussianDistribution sgdistr{};

    MT A(2, 2);
    A << 0.25, 0.75,
         0.75, 3.25;
    Ellipsoid<Point> ell(A);    // origin centered ellipsoid
    GaussianDistribution gdistr(ell);

    // 3. Exponential
    NT variance = 1.0;
    auto c = GetDirection<Point>::apply(HP.dimension(), rng, false);
    ExponentialDistribution edistr(c, variance);

    // 4. LogConcave

    std::pair<Point, NT> inner_ball = HP.ComputeInnerBall();
    Point x0 = inner_ball.first;

    // Reflective HMC and Remmannian HMC are using slightly different functor interfaces
    // TODO: check if this could be unified

    using NegativeGradientFunctorR = CustomFunctor::GradientFunctor<Point>;
    using NegativeLogprobFunctorR = CustomFunctor::FunctionFunctor<Point>;
    CustomFunctor::parameters<NT, Point> params_r(x0);
    NegativeGradientFunctorR gr(params_r);
    NegativeLogprobFunctorR fr(params_r);
    LogConcaveDistribution logconcave_reflective(gr, fr, params_r.L);

    using NegativeGradientFunctor = GaussianFunctor::GradientFunctor<Point>;
    using NegativeLogprobFunctor = GaussianFunctor::FunctionFunctor<Point>;
    using HessianFunctor = GaussianFunctor::HessianFunctor<Point>;
    GaussianFunctor::parameters<NT, Point> params(x0, 0.5, 1);
    NegativeGradientFunctor g(params);
    NegativeLogprobFunctor f(params);
    HessianFunctor h(params);
    LogConcaveDistribution logconcave_crhmc(g, f, h, params.L);

    LogConcaveDistribution logconcave_ref_gaus(g, f, params.L);

    // Sampling

    using NT = double;
    using MT = Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>;
    using VT = Eigen::Matrix<NT,Eigen::Dynamic,1>;

    int rnum = 20;
    int nburns = 5;
    int walk_len = 2;

    MT samples(HP.dimension(), rnum);

    // 1. the eigen matrix interface
    std::cout << "uniform" << std::endl;
    sample_points_eigen_matrix(HP, q, abill_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, abill_walk_custom, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, ball_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, cdhr_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, dikin_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, john_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, rdhr_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, vaidya_walk, udistr, rng, walk_len, rnum, nburns);

    std::cout << "shperical gaussian" << std::endl;
    sample_points_eigen_matrix(HP, q, gball_walk, sgdistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, gcdhr_walk, sgdistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, grdhr_walk, sgdistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, ghmc_walk, sgdistr, rng, walk_len, rnum, nburns);

    std::cout << "general gaussian" << std::endl;
    sample_points_eigen_matrix(HP, q, gbill_walk, gdistr, rng, walk_len, rnum, nburns);

    std::cout << "exponential" << std::endl;
    sample_points_eigen_matrix(HP, q, ehmc_walk, edistr, rng, walk_len, rnum, nburns);

    std::cout << "logconcave" << std::endl;
    sample_points_eigen_matrix(HP, q, hmc_walk, logconcave_reflective, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, nhmc_walk, logconcave_reflective, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, nhmc_walk, logconcave_ref_gaus, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, uld_walk, logconcave_ref_gaus, rng, walk_len, rnum, nburns);

    sample_points_eigen_matrix(HP, q, crhmc_walk, logconcave_crhmc, rng, walk_len, rnum, nburns);
    // The following will compile but segfauls since walk and distribution are not compatible
    //sample_points_eigen_matrix(HP, q, nhmc_walk, logconcave_crhmc, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(problem, q, crhmc_walk, logconcave_crhmc, rng, walk_len, rnum, nburns);


    std::cout << "fix the following" << std::endl;
    // TODO: fix
    // Does not converge because of the starting point
    // Also ess returns rnum instead of 0
    sample_points_eigen_matrix(HP, q, bill_walk, udistr, rng, walk_len, rnum, nburns);

    // Does not compile because of walk-distribution combination
    //sample_points_eigen_matrix(HP, q, abill_walk, gdistr, rng, walk_len, rnum, nburns);

    std::cout << "std::vector interface" << std::endl;
    // 2. the std::vector interface
    std::vector<Point> points;
    sample_points(HP, q, cdhr_walk, udistr, rng, walk_len, rnum, nburns, points);
    for (auto& point : points)
    {
    //    std::cout << point.getCoefficients().transpose() << "\n";
    }

    // 3. the old interface
    // different billiard walks
    typedef BilliardWalk::template Walk<HPolytopeType, RNGType> BilliardWalkType;
    typedef AcceleratedBilliardWalk::template Walk<HPolytopeType, RNGType> AcceleratedBilliardWalkType;
    typedef RandomPointGenerator<AcceleratedBilliardWalkType> Generator;
    std::vector<Point> randPoints;
    PushBackWalkPolicy push_back_policy;
    Generator::apply(HP, q, rnum, walk_len, randPoints, push_back_policy, rng);
    for (auto& point : randPoints)
    {
    //    std::cout << point.getCoefficients().transpose() << "\n";
    }

/*
    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = HP.dimension();
    Point StartingPoint(d);
    std::list<Point> randPoints;

//    gaussian_sampling<AcceleratedBilliardWalk>(points, HP, rng, walkL, numpoints, 1.0,
    //                               StartingPoint, nburns);

    double variance = 1.0;

    Point c(HP.dimension());
    HP.set_InnerBall(std::pair<Point,double>(Point(HP.dimension()), 1.0));
    c = GetDirection<Point>::apply(HP.dimension(), rng, false);
    //ExponentialHamiltonianMonteCarloExactWalk
    exponential_sampling<NutsHamiltonianMonteCarloWalk>(points, HP, rng, walkL, numpoints, c, variance,
                                StartingPoint, nburns);

*/
    return 0;
}