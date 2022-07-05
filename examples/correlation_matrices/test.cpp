#include "matrix_operations/EigenvaluesProblems.h"
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "known_polytope_generators.h"
#include "sampling/sampling.hpp"

#include "diagnostics/univariate_psrf.hpp"

template <typename NT, typename WalkType>
MT call_test(int n, unsigned int num_points=1000, unsigned int walk_len=10, unsigned int nreflex = 10){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    auto start = std::chrono::steady_clock::now();

    MT samples = sample_correlation_matrices<NT, WalkType>(n);

    auto end = std::chrono::steady_clock::now();

    double time = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 2.2);
    return samples;

    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    uniform_sampling<WalkType>(randPoints, P, rng, walkL, numpoints,
                               StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }
}

// Create a TEST_CASE for each random walk

TEST_CASE("hmc") {
    std::cout << "--- Testing HMC" << std::endl;

    int n = 4;
    std::vector<Point> points = call_test();
    for(std::vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
        std::cout << (*it).getCoefficients() << std::endl;
}