#include <iostream>
#include <vector>
#include <fstream>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "sampling/random_point_generators.hpp"
#include "random_walks/random_walks.hpp"
#include "preprocess/max_inscribed_ellipsoid.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "convex_bodies/ellipsoid.h"
#include "convex_bodies/hpolytope.h"


typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef HPolytope <Point> HPOLYTOPE;

// walk policy
PushBackWalkPolicy push_back_policy;

// different billiard walks
typedef BilliardWalk::template Walk<HPOLYTOPE, RNGType> BilliardWalkType;
typedef AcceleratedBilliardWalk::template Walk<HPOLYTOPE, RNGType> AcceleratedBilliardWalkType;
typedef GaussianAcceleratedBilliardWalk::template Walk<HPOLYTOPE, RNGType> GaussianAcceleratedWalkType;


void write_to_file(std::string filename, std::vector<Point> const& randPoints) {
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); //save and redirect
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf); //reset to standard output again
}


void sample_using_uniform_billiard_walk(HPOLYTOPE& HP, RNGType& rng, unsigned int walk_len, unsigned int num_points) {
    std::string filename = "uniform_billiard_walk.txt";
    typedef RandomPointGenerator<BilliardWalkType> Generator;

    std::vector<Point> randPoints;
    Point q(HP.dimension()); // origin
    Generator::apply(HP, q, num_points, walk_len,
                     randPoints, push_back_policy, rng);
    write_to_file(filename, randPoints);
}


void sample_using_uniform_accelerated_billiard_walk(HPOLYTOPE& HP, RNGType& rng, unsigned int walk_len, unsigned int num_points) {
    std::string filename = "uniform_accelerated_billiard_walk.txt";
    typedef RandomPointGenerator<AcceleratedBilliardWalkType> Generator;

    std::vector<Point> randPoints;
    Point q(HP.dimension()); // origin
    Generator::apply(HP, q, num_points, walk_len,
                     randPoints, push_back_policy, rng);
    write_to_file(filename, randPoints);
}


void sample_using_gaussian_billiard_walk(HPOLYTOPE& HP, RNGType& rng, unsigned int walk_len, unsigned int num_points) {
    std::string filename = "gaussian_billiard_walk.txt";
    typedef MultivariateGaussianRandomPointGenerator<GaussianAcceleratedWalkType> Generator;

    std::vector<Point> randPoints;
    Point q(HP.dimension()); // origin

    // ----------- Get inscribed ellipsoid --------------------------------
    unsigned int max_iter = 150;
    NT tol = std::pow(10, -6.0), reg = std::pow(10, -4.0);
    VT x0 = q.getCoefficients();
    VT center;
    bool converged;
    std::tuple<MT, VT, NT> ellipsoid = compute_inscribed_ellipsoid<MT, EllipsoidType::VOLUMETRIC_BARRIER>
    (HP.get_mat(), HP.get_vec(), x0, max_iter, tol, reg);
    
    const MT E = get<0>(ellipsoid);
    // --------------------------------------------------------------------

    Generator::apply(HP, q, E, num_points, walk_len,
                     randPoints, push_back_policy, rng);
    write_to_file(filename, randPoints);
}


HPOLYTOPE visualization_example() {
    MT Ab(5, 3);
    Ab <<  -1,  0, 1,     // x >= -1
            0, -1, 1,     // y >= -1
            1, -1, 8,     // x-y <= 8
           -1,  1, 8,     // x-y >= -8
            1,  1, 50;    // x+y <= 50

    MT A = Ab.block(0, 0, 5, 2);
    VT b = Ab.col(2);
    return HPOLYTOPE(2, A, b);
}


int main(int argc, char const *argv[]) {
    unsigned int walk_len = 10;
    unsigned int num_points = 1000;

    HPOLYTOPE HP = visualization_example();
    HP.ComputeInnerBall();

    RNGType rng(HP.dimension());
    sample_using_uniform_billiard_walk(HP, rng, walk_len, num_points);
    sample_using_uniform_accelerated_billiard_walk(HP, rng, walk_len, num_points);
    sample_using_gaussian_billiard_walk(HP, rng, walk_len, num_points);
    return 0;
}
