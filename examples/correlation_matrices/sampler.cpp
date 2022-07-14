#ifndef EIGCORRELATION
    #define EIGCORRELATION
#endif

#include <vector>
#include <chrono>
#include <iostream>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "sampling/random_point_generators.hpp"
#include "random_walks/random_walks.hpp"
#include "convex_bodies/correlation_matrices/corre_spectra.hpp"
#include "diagnostics/univariate_psrf.hpp"

#include "matrix_operations/EigenvaluesProblems.h"
#include "generators/boost_random_number_generator.hpp"
#include "misc.h"
#include "random_walks/random_walks.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "direct_sampler.hpp"

template<typename PointType,
    typename PointList,
    typename RandomNumberGenerator, 
    typename WalkTypePolicy>
void uniform_correlation_sampling(PointList &randPoints,
                   CorreSpectra<PointType> &P,
                   RandomNumberGenerator &rng,
                   const unsigned int &walkL,
                   const unsigned int &num_points,
                   const PointType &starting_point,
                   unsigned int const& nburns){
    typedef typename WalkTypePolicy::template Walk <CorreSpectra<PointType>,
                                                    RandomNumberGenerator> walk;
    PushBackWalkPolicy push_back_policy;
    typedef RandomPointGenerator<walk> RandomPointGenerator;
    
    PointType p = starting_point;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walkL, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, num_points, walkL, randPoints,
                                push_back_policy, rng);
}

template <typename NT, typename WalkType, typename RNGType>
Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> call_test(unsigned int n, unsigned int const num_points, unsigned int walkL, unsigned int nreflex){

    std::cout << "Improved implementation : " << std::endl;

    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef std::vector<Point> PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 
    
    typedef CorreSpectra<Point>     CorreSpectraType;
    typedef RandomPointGenerator<WalkType> Generator;

    std::vector<Point> randPoints;
    CorreSpectraType P(n);

    const unsigned int nburns = 0, d = P.dimension();
    Point startingPoint(d);
    RNGType rng(d);
    
    auto start = std::chrono::steady_clock::now();

    uniform_correlation_sampling<Point, PointList, RNGType, WalkType>(randPoints, P, rng, walkL, num_points, startingPoint, nburns);

    auto end = std::chrono::steady_clock::now();

    double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Improved time : " << time << std::endl;
    // write_to_file<Point>("uniform_billiard_walk.txt", randPoints, time);
    
    MT samples(d,num_points);
    // int j = 0;
    // for (typename std::vector<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); ++rpit, ++j)
    //     samples.col(j) = (*rpit).getCoefficients();

    // VT score = univariate_psrf<NT, VT, MT>(samples);
    // std::cout << "psrf = " << score.maxCoeff() << std::endl;

    // CHECK(score.maxCoeff() < 2.2);
    
    return samples;
}

int main(int argc, char const *argv[]) {
    // srand((unsigned) time(NULL));
    srand(19031999);
    typedef double NT;
    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    
    // BilliardWalk, AcceleratedBilliardWalk, GaussianAcceleratedBilliardWalk
    // GaussianHamiltonianMonteCarloExactWalk
    unsigned int n, num_points = 1000, walkL = 10, nreflex = 10;
    std::cout << "n = ";
    std::cin >> n;
    // MT M = call_test<NT, BilliardWalk, RNGType>(n, num_points, walkL, nreflex);

    // Test direct implementation:

    naive_test<NT, BilliardWalk, RNGType>(n, num_points, walkL, nreflex);

    // naive_test<NT, GaussianHamiltonianMonteCarloExactWalk, RNGType>(n, num_points, walkL, nreflex);

    // std::cout << M << std::endl;
    
    return 0;
}