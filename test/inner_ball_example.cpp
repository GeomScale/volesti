#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"

using NT = double;
using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
using Kernel = Cartesian<NT>;
using Point = typename Kernel::Point;
using Polytope = HPolytope<Point>;

int main() {
    MT A(4,2);
    VT b(4);

    A <<  1,  0,
         -1,  0,
          0,  1,
          0, -1;
    b << 1, 1, 1, 1;   // unit square

    Polytope P(2, A, b);

    auto inner_ball = P.ComputeInnerBall();
    Point c = inner_ball.first;
    NT r = inner_ball.second;

    std::cout << "radius = " << r << "\ncenter = ";
    for (unsigned i = 0; i < P.dimension(); ++i) std::cout << c[i] << " ";
    std::cout << "\n";
}