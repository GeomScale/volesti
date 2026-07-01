// VolEsti (volume computation and sampling library)

// Copyright (c) 2024-2025 Vissarion Fisikopoulos
// Copyright (c) 2024-2025 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <chrono>
#include <boost/random.hpp>
#include "generators/boost_random_number_generator.hpp"
#include "generators/known_polytope_generators.h"
#include "generators/z_polytopes_generators.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/convex_body.h"
#include "generators/convex_bodies_generator.h"
#include "cartesian_geom/cartesian_kernel.h"

template <typename NT>
void call_test_known_generators()
{
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;

    // 1. generate_cube<HPolytope>(dim, false): cube [-1,1]^dim
    {
        unsigned int dim = 2;
        Hpolytope P = generate_cube<Hpolytope>(dim, false);
        CHECK(P.dimension() == dim);
        CHECK(P.num_of_hyperplanes() == 2 * (int)dim);

        Point origin(dim);
        CHECK(P.is_in(origin) == -1);

        Point outside(dim);
        outside.set_coord(0, 1.5);
        outside.set_coord(1, 0.0);
        CHECK(P.is_in(outside) == 0);

        std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
        CHECK(inner_ball.second > NT(0));
    }

    // 2. generate_cube<HPolytope>(dim, false, 2.0): cube [-2,2]^dim
    {
        unsigned int dim = 2;
        Hpolytope P = generate_cube<Hpolytope>(dim, false, 2.0);
        CHECK(P.dimension() == dim);
        CHECK(P.num_of_hyperplanes() == 2 * (int)dim);

        Point origin(dim);
        CHECK(P.is_in(origin) == -1);

        // (1.5, 0) should be inside when scale=2
        Point inside(dim);
        inside.set_coord(0, 1.5);
        CHECK(P.is_in(inside) == -1);

        // (2.5, 0) should be outside
        Point outside(dim);
        outside.set_coord(0, 2.5);
        CHECK(P.is_in(outside) == 0);
    }

    // 3. generate_cube<VPolytope>(dim, true): V-representation cube
    {
        unsigned int dim = 2;
        Vpolytope P = generate_cube<Vpolytope>(dim, true);
        CHECK(P.dimension() == dim);
        CHECK(P.num_of_vertices() == (1 << dim)); // 2^dim = 4
        CHECK(P.num_of_vertices() == 4);

        Point origin(dim);
        CHECK(P.is_in(origin) == -1);
    }

    // 4. generate_cross<HPolytope>(dim, false): L1 ball |x|_1 <= 1
    {
        unsigned int dim = 2;
        Hpolytope P = generate_cross<Hpolytope>(dim, false);
        CHECK(P.dimension() == dim);

        Point origin(dim);
        CHECK(P.is_in(origin) == -1);

        // A basis vector e1 is on the boundary since |e1|_1 = 1 <= 1
        Point basis(dim);
        basis.set_coord(0, 1.0);
        CHECK(P.is_in(basis) == -1);

        // (0.6, 0.6) gives L1-norm 1.2 > 1 so it is outside
        Point outside(dim);
        outside.set_coord(0, 0.6);
        outside.set_coord(1, 0.6);
        CHECK(P.is_in(outside) == 0);
    }

    // 5. generate_simplex<HPolytope>(dim, false):
    //    Constraints: x_i <= 0, -sum(x_i) <= 1  i.e. sum(x_i) >= -1
    //    Vertices: (0,0), (-1,0), (0,-1) for dim=2
    {
        unsigned int dim = 2;
        Hpolytope P = generate_simplex<Hpolytope>(dim, false);
        CHECK(P.dimension() == dim);

        Point origin(dim);
        CHECK(P.is_in(origin) == -1);

        // (1, 0) violates x1 <= 0 and is outside.
        Point outside_pos(dim);
        outside_pos.set_coord(0, 1.0);
        CHECK(P.is_in(outside_pos) == 0);

        // (-0.5, -0.5) satisfies x_i <= 0 and sum >= -1, so it's inside
        Point inside_neg(dim);
        inside_neg.set_coord(0, -0.5);
        inside_neg.set_coord(1, -0.5);
        CHECK(P.is_in(inside_neg) == -1);
    }

    // 6. generate_prod_simplex<HPolytope>(dim):
    //    Product of two simplices, dimension = 2*dim
    {
        unsigned int dim = 2;
        Hpolytope P = generate_prod_simplex<Hpolytope>(dim);
        CHECK(P.dimension() == 2 * dim);

        Point origin(2 * dim);
        CHECK(P.is_in(origin) == -1);

        // Vpoly=true is not supported; returns a default-constructed polytope with 0 hyperplanes
        Hpolytope Perr = generate_prod_simplex<Hpolytope>(dim, true);
        CHECK(Perr.num_of_hyperplanes() == 0);
    }

    // 7. generate_skinny_cube<HPolytope>(dim):
    //    Cube with first axis scaled to [-100, 100], others at [-1, 1]
    {
        unsigned int dim = 3;
        Hpolytope P = generate_skinny_cube<Hpolytope>(dim);
        CHECK(P.dimension() == dim);

        Point origin(dim);
        CHECK(P.is_in(origin) == -1);

        // (0, 2, 0) should be outside because x2 in [-1, 1]
        Point outside(dim);
        outside.set_coord(1, 2.0);
        CHECK(P.is_in(outside) == 0);
    }

    // 8. generate_birkhoff<HPolytope>(n):
    //    Birkhoff polytope of doubly stochastic n x n matrices
    //    dim = n^2 - 2n + 1, num_hyperplanes = n^2
    {
        unsigned int n = 3;
        Hpolytope P = generate_birkhoff<Hpolytope>(n);
        unsigned int expected_dim = n * n - 2 * n + 1; // 9 - 6 + 1 = 4
        CHECK(P.dimension() == expected_dim);
        CHECK(P.num_of_hyperplanes() == (int)(n * n)); // 9
    }
}

TEST_CASE("known_polytope_generators")
{
    call_test_known_generators<double>();
}
