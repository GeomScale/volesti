#include "doctest.h"

#include <Eigen/Eigen>

#include "convex_bodies/simplexintersectball.h"
#include "cartesian_geom/cartesian_kernel.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef SimplexIntersectBall<Point> SimplexBall;

TEST_CASE("simplexintersectball_basic_membership")
{
    unsigned int d = 2;

    // Simplex in R^2:
    // x >= 0, y >= 0, x + y <= 1
    MT A(3, 2);
    A << -1,  0,
          0, -1,
          1,  1;

    VT b(3);
    b << 0, 0, 1;

    // Vertices stored column-wise: (0,0), (1,0), (0,1)
    MT V(2, 3);
    V << 0, 1, 0,
         0, 0, 1;

    VT x0(2);
    x0 << 0, 0;

    SimplexBall K(d, A, b, V, x0);

    CHECK(K.dimension() == 2);
    CHECK(K.num_of_hyperplanes() == 3);

    VT inside_vec(2);
    inside_vec << 0.2, 0.2;
    Point inside(inside_vec);

    VT outside_simplex_vec(2);
    outside_simplex_vec << 0.8, 0.8;
    Point outside_simplex(outside_simplex_vec);

    VT outside_ball_vec(2);
    outside_ball_vec << 1.2, 0.0;
    Point outside_ball(outside_ball_vec);

    CHECK(K.is_in(inside) == -1);
    CHECK(K.is_in(outside_simplex) == 0);
    CHECK(K.is_in(outside_ball) == 0);
}

TEST_CASE("simplexintersectball_line_intersection")
{
    unsigned int d = 2;

    MT A(3, 2);
    A << -1,  0,
          0, -1,
          1,  1;

    VT b(3);
    b << 0, 0, 1;

    MT V(2, 3);
    V << 0, 1, 0,
         0, 0, 1;

    VT x0(2);
    x0 << 0, 0;

    SimplexBall K(d, A, b, V, x0);

    VT r_vec(2);
    r_vec << 0.2, 0.2;
    Point r(r_vec);

    VT v_vec(2);
    v_vec << 1.0, 0.0;
    Point v(v_vec);

    std::pair<NT, NT> interval = K.line_intersect(r, v);

    CHECK(interval.first == doctest::Approx(0.6));
    CHECK(interval.second == doctest::Approx(-0.2));
}