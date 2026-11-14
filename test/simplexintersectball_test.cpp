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

SimplexBall make_3d_tetrahedron_unit_ball()
{
    unsigned int d = 3;

    // Tetrahedron in R^3:
    // x >= 0, y >= 0, z >= 0, x + y + z <= 2
    MT A(4, 3);
    A << -1, 0, 0,
        0, -1, 0,
        0, 0, -1,
        1, 1, 1;

    VT b(4);
    b << 0, 0, 0, 2;

    // Vertices stored column-wise:
    // (0,0,0), (2,0,0), (0,2,0), (0,0,2)
    MT V(3, 4);
    V << 0, 2, 0, 0,
        0, 0, 2, 0,
        0, 0, 0, 2;

    VT x0(3);
    x0 << 0, 0, 0;

    return SimplexBall(d, A, b, V, x0);
}

TEST_CASE("simplexintersectball_basic_membership")
{
    unsigned int d = 2;

    // Simplex in R^2:
    // x >= 0, y >= 0, x + y <= 1
    MT A(3, 2);
    A << -1, 0,
        0, -1,
        1, 1;

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
    A << -1, 0,
        0, -1,
        1, 1;

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

TEST_CASE("simplexintersectball_3d_tetrahedron_membership")
{
    SimplexBall K = make_3d_tetrahedron_unit_ball();

    CHECK(K.dimension() == 3);
    CHECK(K.num_of_hyperplanes() == 4);

    VT inside_vec(3);
    inside_vec << 0.3, 0.3, 0.3;
    Point inside(inside_vec);

    VT outside_simplex_vec(3);
    outside_simplex_vec << -0.1, 0.2, 0.2;
    Point outside_simplex(outside_simplex_vec);

    VT outside_ball_vec(3);
    outside_ball_vec << 1.2, 0.0, 0.0;
    Point outside_ball(outside_ball_vec);

    CHECK(K.is_in(inside) == -1);
    CHECK(K.is_in(outside_simplex) == 0);
    CHECK(K.is_in(outside_ball) == 0);
}

TEST_CASE("simplexintersectball_3d_tetrahedron_gc_intersection")
{
    SimplexBall K = make_3d_tetrahedron_unit_ball();

    NT inv_sqrt3 = NT(1) / std::sqrt(NT(3));
    NT inv_sqrt2 = NT(1) / std::sqrt(NT(2));

    VT r_vec(3);
    r_vec << inv_sqrt3, inv_sqrt3, inv_sqrt3;
    Point r(r_vec);

    VT v_vec(3);
    v_vec << inv_sqrt2, -inv_sqrt2, 0;
    Point v(v_vec);

    VT Ar;
    VT Av;

    std::pair<NT, NT> interval = K.gc_intersect(r, v, Ar, Av);

    NT alpha = std::atan(std::sqrt(NT(2) / NT(3)));

    CHECK(interval.first == doctest::Approx(alpha));
    CHECK(interval.second == doctest::Approx(-alpha));
}

TEST_CASE("simplexintersectball_3d_tetrahedron_gc_intersection_positive")
{
    SimplexBall K = make_3d_tetrahedron_unit_ball();

    NT inv_sqrt3 = NT(1) / std::sqrt(NT(3));
    NT inv_sqrt2 = NT(1) / std::sqrt(NT(2));

    VT r_vec(3);
    r_vec << inv_sqrt3, inv_sqrt3, inv_sqrt3;
    Point r(r_vec);

    VT v_vec(3);
    v_vec << inv_sqrt2, -inv_sqrt2, 0;
    Point v(v_vec);

    VT Ar;
    VT Av;

    std::pair<NT, int> hit = K.gc_intersect_positive(r, v, Ar, Av);

    NT alpha = std::atan(std::sqrt(NT(2) / NT(3)));

    CHECK(hit.first == doctest::Approx(alpha));
    CHECK(hit.second == 1);
}

TEST_CASE("simplexintersectball_incremental_rotation_matches_fresh")
{
    SimplexBall K = make_3d_tetrahedron_unit_ball();

    NT inv_sqrt3 = NT(1) / std::sqrt(NT(3));
    NT inv_sqrt2 = NT(1) / std::sqrt(NT(2));

    VT r0(3);
    r0 << inv_sqrt3, inv_sqrt3, inv_sqrt3;
    VT v0(3);
    v0 << inv_sqrt2, -inv_sqrt2, 0;
    Point pr0(r0), pv0(v0);

    VT Ar, Av;
    K.gc_intersect(pr0, pv0, Ar, Av);

    NT lambda = 0.3;
    std::pair<NT, NT> inc = K.gc_intersect(pr0, pv0, Ar, Av, lambda);

    VT r1 = std::cos(lambda) * r0 + std::sin(lambda) * v0;
    Point pr1(r1);
    VT Ar2, Av2;
    std::pair<NT, NT> fresh = K.gc_intersect(pr1, pv0, Ar2, Av2);

    CHECK(inc.first == doctest::Approx(fresh.first));
    CHECK(inc.second == doctest::Approx(fresh.second));
}

TEST_CASE("simplexintersectball_3d_tetrahedron_reflection")
{
    SimplexBall K = make_3d_tetrahedron_unit_ball();

    NT inv_sqrt2 = NT(1) / std::sqrt(NT(2));
    // Point on the unit sphere and on the facet y = 0
    VT p(3);
    p << inv_sqrt2, 0, inv_sqrt2;

    // Tangent direction pointing outside through y < 0
    VT v(3);
    v << 0, -1, 0;

    MT projector = MT::Identity(3, 3) - p * p.transpose();

    int facet = 1; // y = 0 facet, represented by -y <= 0

    K.compute_reflection(v, p, projector, facet);

    CHECK(v(0) == doctest::Approx(0.0));
    CHECK(v(1) == doctest::Approx(1.0));
    CHECK(v(2) == doctest::Approx(0.0));
}

TEST_CASE("simplexintersectball_3d_tetrahedron_gc_all_roots")
{
    SimplexBall K = make_3d_tetrahedron_unit_ball();

    NT inv_sqrt3 = NT(1) / std::sqrt(NT(3));
    NT inv_sqrt2 = NT(1) / std::sqrt(NT(2));

    VT r_vec(3);
    r_vec << inv_sqrt3, inv_sqrt3, inv_sqrt3;
    Point r(r_vec);

    VT v_vec(3);
    v_vec << inv_sqrt2, -inv_sqrt2, 0;
    Point v(v_vec);

    VT Ar;
    VT Av;

    std::pair<VT, VT> roots = K.gc_intersect_all_roots(r, v, Ar, Av);

    NT alpha = std::atan(std::sqrt(NT(2) / NT(3)));

    CHECK(roots.first.rows() >= 1);
    CHECK(roots.second.rows() >= 1);

    bool found_negative_alpha = false;
    for (int i = 0; i < roots.first.rows(); ++i)
    {
        if (std::abs(roots.first(i) + alpha) < NT(1e-08))
        {
            found_negative_alpha = true;
        }
    }

    bool found_positive_alpha = false;
    for (int i = 0; i < roots.second.rows(); ++i)
    {
        if (std::abs(roots.second(i) - alpha) < NT(1e-08))
        {
            found_positive_alpha = true;
        }
    }

    CHECK(found_negative_alpha);
    CHECK(found_positive_alpha);
}