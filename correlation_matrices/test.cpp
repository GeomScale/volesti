#include "doctest.h"
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
void call_test(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    unsigned int n = 10;

    auto start = std::chrono::steady_clock::now();

    MT samples = sample_correlation_matrices<NT, WalkType>(n);

    auto end = std::chrono::steady_clock::now();

    double time = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 2.2);
}

// Create a TEST_CASE for each random walk

TEST_CASE("hmc") {
    std::cout << "--- Testing HMC" << std::endl;

    int n = 4;
    std::vector<Point> points = call_test();
    for(std::vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
        std::cout << (*it).getCoefficients() << std::endl;
}