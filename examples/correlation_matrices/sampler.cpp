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

#include "uniform_correlation_matrices.hpp"
#include "direct_sampler.hpp"
#include "test.hpp"

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
    unsigned int n, num_points = 100, walkL = 10;
    std::cout << "n = ";
    std::cin >> n;
    
    new_test<NT, BilliardWalk, RNGType>(n, num_points, walkL);

    naive_uniform_test<NT, BilliardWalk, RNGType>(n, num_points, walkL);

    // naive_test<NT, GaussianHamiltonianMonteCarloExactWalk, RNGType>(n, num_points, walkL, nreflex);

    // std::cout << M << std::endl;
    
    return 0;
}