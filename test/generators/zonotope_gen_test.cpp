// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <chrono>
#include <boost/random.hpp>

#include "generators/z_polytopes_generators.h"
#include "convex_bodies/zpolytope.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"

template <typename NT>
void call_test_zono_gens() {

    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef Zonotope<Point> Zono;
    typedef boost::mt19937 RNGType;

    int dim = 2;
    int m = 5;
    double seed = 42.0;
    Point origin(dim);

    // Test gen_zonotope_gaussian
    Zono P1 = gen_zonotope_gaussian<Zono, RNGType>(dim, m, seed);
    CHECK(P1.dimension() == static_cast<unsigned int>(dim));
    CHECK(P1.num_of_generators() == m);
    CHECK(P1.is_in(origin) == -1);

    // Reproducibility: same seed should give identical zonotope
    Zono P1_copy = gen_zonotope_gaussian<Zono, RNGType>(dim, m, seed);
    CHECK(P1_copy.dimension() == static_cast<unsigned int>(dim));
    CHECK(P1_copy.num_of_generators() == m);
    CHECK(P1_copy.is_in(origin) == -1);

    // Test gen_zonotope_uniform
    Zono P2 = gen_zonotope_uniform<Zono, RNGType>(dim, m, seed);
    CHECK(P2.dimension() == static_cast<unsigned int>(dim));
    CHECK(P2.num_of_generators() == m);
    CHECK(P2.is_in(origin) == -1);

    // Reproducibility
    Zono P2_copy = gen_zonotope_uniform<Zono, RNGType>(dim, m, seed);
    CHECK(P2_copy.dimension() == static_cast<unsigned int>(dim));
    CHECK(P2_copy.num_of_generators() == m);
    CHECK(P2_copy.is_in(origin) == -1);

    // Test gen_zonotope_exponential
    Zono P3 = gen_zonotope_exponential<Zono, RNGType>(dim, m, seed);
    CHECK(P3.dimension() == static_cast<unsigned int>(dim));
    CHECK(P3.num_of_generators() == m);
    CHECK(P3.is_in(origin) == -1);

    // Reproducibility
    Zono P3_copy = gen_zonotope_exponential<Zono, RNGType>(dim, m, seed);
    CHECK(P3_copy.dimension() == static_cast<unsigned int>(dim));
    CHECK(P3_copy.num_of_generators() == m);
    CHECK(P3_copy.is_in(origin) == -1);
}

TEST_CASE("zonotope_generators") {
    call_test_zono_gens<double>();
}
