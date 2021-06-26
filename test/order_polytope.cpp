// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021- Vaibhav Thakkar

// Contributed by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include "cartesian_geom/point.h"
#include "poset.h"
#include "orderpolytope.h"
#include "misc.h"
#include "random.hpp"



// template <typename NT>
// void test_values(NT computed, NT expected)
// {
//     std::cout << "Computed value " << computed << std::endl;
//     std::cout << "Expected value = " << expected << std::endl;
//     std::cout << "Relative error (expected) = "
//               << std::abs((computed-expected)/expected) << std::endl;
//     CHECK(std::abs((volume - expected)/expected) < 0.00001);
// }


// template <class Polytope>
// void test_volume(Polytope &HP,
//                  double const& expectedBall,
//                  double const& expectedCDHR,
//                  double const& expectedRDHR,
//                  double const& expectedBilliard,
//                  double const& exact,
//                  bool birk = false)
// {
//     typedef typename Polytope::PointType Point;
//     typedef typename Point::FT NT;

//     // Setup the parameters
//     int walk_len = 10 + HP.dimension()/10;
//     NT e=0.1, volume;

//     // Estimate the volume
//     std::cout << "Number type: " << typeid(NT).name() << std::endl;
//     typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

//     //TODO: low accuracy in high dimensions
//     if (!birk) {
//         volume = volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len).second;
//         test_values(volume, expectedBall, exact);
//     }

//     volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, e, walk_len).second;
//     test_values(volume, expectedCDHR, exact);

//     volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, e, walk_len).second;
//     test_values(volume, expectedRDHR, exact);

//     volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
//     test_values(volume, expectedBilliard, exact);
// }


// template <typename NT>
// void call_test_reflection()
// {
//     typedef Cartesian<NT>    Kernel;
//     typedef typename Kernel::Point    Point;
//     typedef HPolytope<Point> Hpolytope;
//     Hpolytope P;

//     typedef BoostRandomNumberGenerator<boost::mt19937, NT, 123> RNGType;

//     std::cout << "--- Testing volume of H-birk3" << std::endl;
//     P = generate_birkhoff<Hpolytope>(3);
//     test_volume(P, 0.114343, 0.125548, 0.113241, 0.112446, 0.125, true);

//     std::cout << "--- Testing volume of H-birk4" << std::endl;
//     P = generate_birkhoff<Hpolytope>(4);
//     test_volume(P, 0.00112956, 0.00109593, 0.00108152, 0.000845192,
//                 0.000970018, true);

//     std::cout << "--- Testing volume of H-birk5" << std::endl;
//     P = generate_birkhoff<Hpolytope>(5);
//     test_volume(P,
//                 1.97968e-07,
//                 1.73729e-07,
//                 1.39042e-07,
//                 3.24308e-07,
//                 0.000000225, 
//                 true);

//     std::cout << "--- Testing volume of H-birk6" << std::endl;
//     P = generate_birkhoff<Hpolytope>(6);
//     test_volume(P,
//                 7.84351e-13,
//                 6.10783e-13,
//                 5.05917e-13,
//                 6.62349e-13,
//                 9.455459196 * std::pow(10,-13), 
//                 true);
// }

// template <typename NT>
// void call_test_line_intersect() {
//     typedef Cartesian<NT>    Kernel;
//     typedef typename Kernel::Point    Point;

//     typedef HPolytope<Point> Hpolytope;
//     Hpolytope P;

//     std::cout << "--- Testing volume of H-prod_simplex5" << std::endl;
//     P = generate_prod_simplex<Hpolytope>(5);
//     test_volume(P,
//                 6.40072 * std::pow(10,-5),
//                 6.69062 * std::pow(10,-5),
//                 6.20744e-05,
//                 6.31986 * std::pow(10,-5),
//                 std::pow(1.0 / factorial(5.0), 2));

//     std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
//     P = generate_prod_simplex<Hpolytope>(10);
//     test_volume(P,
//                 6.83631 * std::pow(10,-14),
//                 8.19581 * std::pow(10,-14),
//                 9.35005e-14,
//                 6.57309e-14,
//                 std::pow(1.0 / factorial(10.0), 2));

//     std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
//     P = generate_prod_simplex<Hpolytope>(15);
//     test_volume(P,
//                 3.85153e-25,
//                 9.33162 * std::pow(10,-25),
//                 3.95891e-25,
//                 5.72542e-25,
//                 std::pow(1.0 / factorial(15.0), 2));
// }

// template <typename NT>
// void call_test_is_in() {
//     typedef Cartesian<NT>    Kernel;
//     typedef typename Kernel::Point    Point;

//     typedef HPolytope<Point> Hpolytope;
//     Hpolytope P;

//     std::cout << "--- Testing volume of H-simplex10" << std::endl;
//     P = generate_simplex<Hpolytope>(10, false);
//     test_volume(P,
//                 3.90133e-07,
//                 2.90617 * std::pow(10,-7),
//                 2.93392 * std::pow(10,-7),
//                 3.00286e-07,
//                 1.0 / factorial(10.0));

//     std::cout << "--- Testing volume of H-simplex20" << std::endl;
//     P = generate_simplex<Hpolytope>(20, false);
//     test_volume(P,
//                 6.52535e-19,
//                 4.14182 * std::pow(10,-19),
//                 4.5877e-19,
//                 4.54245e-19,
//                 1.0 / factorial(20.0));

//     std::cout << "--- Testing volume of H-simplex30" << std::endl;
//     P = generate_simplex<Hpolytope>(30, false);
//     test_volume(P,
//                 2.5776 * std::pow(10,-33),
//                 3.5157 * std::pow(10,-33),
//                 2.74483e-33,
//                 3.08769e-33,
//                 1.0 / factorial(30.0));
// }

template <typename NT>
void call_test_basics() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point Point;
    typedef typename OrderPolytope<Point>::VT VT;
    typedef typename Poset::RT RT;
    typedef typename Poset::RV RV;


    // Create Poset, 4 elements, a0 <= a1, a0 <= a2, a1 <= a3
    RV poset_data{{0, 1}, {0, 2}, {1, 3}};
    Poset poset(4, poset_data);
    CHECK(poset.num_elem() == 4);
    CHECK(poset.num_relations() == 3);

    
    // Initialize order polytope from the poset
    OrderPolytope<Point> OP(poset);
    unsigned int d = OP.dimension(), m = OP.num_hyperplanes();
    CHECK(d == 4);
    CHECK(m == 2*4 + 3);


    VT expected_dist_vector = Eigen::MatrixXd::Zero(m, 1);
    expected_dist_vector.block(d, 0, d, 1) = Eigen::MatrixXd::Constant(d, 1, 1.0);
    VT ret_dists_vector = Eigen::Map< VT >(OP.get_dists(0.0).data(), m);
    CHECK( (expected_dist_vector - ret_dists_vector).norm() == 0 );
    
    CHECK(OP.is_in(Point(4, {0.0, 0.5, 1.0, 1.0})) == -1);
    CHECK(OP.is_in(Point(4, {1.0, 0.5, 1.0, 1.0})) == 0);   // a0 <= a1 violated
    CHECK(OP.is_in(Point(4, {0.5, 0.5, 0.0, 1.0})) == 0);   // a0 <= a2 violated
    CHECK(OP.is_in(Point(4, {-0.1, 0.5, 1.0, 1.0})) == 0);  // a0 >= 0 violated
    CHECK(OP.is_in(Point(4, {1.0, 0.5, 1.0, 1.1})) == 0);   // a3 <= 1 violated
}


TEST_CASE("basics") {
    call_test_basics<double>();
}

TEST_CASE("line_intersect") {
    // call_test_line_intersect<double>();
}

TEST_CASE("reflection") {
    // call_test_reflection<double>();
}

TEST_CASE("vec_mult") {
    // call_test_vec_mult<double>();
}
