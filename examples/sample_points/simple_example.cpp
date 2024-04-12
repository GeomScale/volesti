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

template <typename Walk, typename Distribution>
void sample_points_eigen_matrix(HPolytopeType const& HP, Point const& q, Walk const& walk,
                                Distribution const& distr, RNGType rng, int walk_len, int rnum,
                                int nburns)
{
    using NT = double;
    using MT = Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>;
    using VT = Eigen::Matrix<NT,Eigen::Dynamic,1>;

    MT samples(HP.dimension(), rnum);

    sample_points(HP, q, walk, distr, rng, walk_len, rnum, nburns, samples);

    // sample stats
    unsigned int min_ess;
    auto score = effective_sample_size<NT, VT, MT>(samples, min_ess);
    std::cout << "ess=" << min_ess << std::endl;
    //print_diagnostics<NT, VT, MT>(samples, min_ess, std::cerr);
}

int main() {
    // Generating a 3-dimensional cube centered at origin
    HPolytopeType HP = generate_cube<HPolytopeType>(10, false);
    std::cout<<"Polytope: \n";
    HP.ComputeInnerBall();
    //HP.print();
    //std::cout<<"\n";

    // Setup parameters for sampling
    Point q(HP.dimension());
    RNGType rng(HP.dimension());

    // NEW INTERFACE Sampling
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

    UniformDistribution udistr{};
    GaussianDistribution gdistr{};

    using NT = double;
    using MT = Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>;
    using VT = Eigen::Matrix<NT,Eigen::Dynamic,1>;

    int rnum = 100;
    int nburns = 50;
    int walk_len = 10;

    MT samples(HP.dimension(), rnum);

    // 1. the eigen matrix interface
    sample_points_eigen_matrix(HP, q, abill_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, abill_walk_custom, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, ball_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, cdhr_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, dikin_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, john_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, rdhr_walk, udistr, rng, walk_len, rnum, nburns);
    sample_points_eigen_matrix(HP, q, vaidya_walk, udistr, rng, walk_len, rnum, nburns);

    // TODO: fix
    // Does not converge because of the starting point
    // Also ess returns rnum instead of 0
    sample_points_eigen_matrix(HP, q, bill_walk, udistr, rng, walk_len, rnum, nburns);

    // Does not compile because of walk-distribution combination
    //sample_points_eigen_matrix(HP, q, abill_walk, gdistr, rng, walk_len, rnum, nburns);

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