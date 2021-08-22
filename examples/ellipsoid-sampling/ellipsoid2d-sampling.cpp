#include <iostream>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "ellipsoid.h"
#include "sampling/ellipsoid.hpp"
#include "random_walks/random_walks.hpp"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

int main(int argc, char const *argv[]) {
    unsigned int dim = 2;
    MT A(2, 2);
    A << 0.25, 0.75,
         0.75, 3.25;

    Ellipsoid<Point> ell(A);    // origin centered ellipsoid
    int num_points = 1000;
    Point p(dim);
    RNGType rng(dim);

    for (int i=0; i<num_points; ++i) {
        p = GetPointInDellipsoid<Point>::apply(dim, ell, rng);
        p.print();
    }
    return 0;
}

