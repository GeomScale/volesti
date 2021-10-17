// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "convex_bodies/hpolytope.h"
#include "Eigen/Eigen"
#include "doctest.h"
#include "LV_MC_integration.hpp"
#include "ode_solvers/oracle_functors.hpp"
#include "generators/known_polytope_generators.h"
#include "boost_random_number_generator.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

template <typename NT>
void test_values (NT computed, NT expected, NT exact) {
	std::cout << "Computed integration value = " << computed << std::endl;
	std::cout << "Expected integration value = " << expected << std::endl;
	std::cout << "Exact integration value = " << exact << std::endl;
	std::cout << "Relative error (expected) = " << std::abs((computed - expected)/expected) << std::endl;
	std::cout << "Relative error (exact) = " << std::abs((computed - exact)/exact) << std::endl ;
	CHECK(((std::abs((computed - expected)/expected) < 0.00001) || (std::abs((computed - exact)/exact) < 0.2)));
}

template <typename NT>
void cubes_test_lovasz_vempala_integrate (const int dimension, NT expected, NT exact, NT beta = 1.0, volumetype voltype = SOB) {

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef IsotropicQuadraticFunctor::parameters<NT> IsoParameters;

	IsoParameters params;
	params.alpha = (NT)2;
	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);
	
	HPOLYTOPE HP;
	HP = generate_cube <HPOLYTOPE> (dimension, false);
	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();;
	Point x0 = inner_ball.first;
	NT integral_value = lovasz_vempala_integrate 
	  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, HP, x0, beta, voltype, 5, 0.1);
	
	test_values(integral_value, expected, exact);
}

template <typename NT>
void simplices_test_lovasz_vempala_integrate (const int dimension, NT expected, NT exact, NT beta = 1.0, volumetype voltype = SOB) {

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef IsotropicQuadraticFunctor::parameters<NT> IsoParameters;

	IsoParameters params;
	params.alpha = (NT)2;
	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	HPOLYTOPE HP;
	HP = generate_simplex <HPOLYTOPE> (dimension, false);
	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;
	NT integral_value = lovasz_vempala_integrate 
	  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, HP, x0, beta, voltype, 5, 0.1);
	
	test_values(integral_value, expected, exact);
}

template <typename NT>
void prod_simplices_test_lovasz_vempala_integrate (const int dimension, NT expected, NT exact, NT beta = 1.0, volumetype voltype = SOB) {

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef IsotropicQuadraticFunctor::parameters<NT> IsoParameters;

	IsoParameters params;
	params.alpha = (NT)2;
	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);
	HPOLYTOPE HP;
	HP = generate_prod_simplex <HPOLYTOPE> (dimension, false);
	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;
	NT integral_value = lovasz_vempala_integrate 
	  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, HP, x0, beta, voltype, 5, 0.1);
	
	test_values(integral_value, expected, exact);
}

template <typename NT>
void cross_test_lovasz_vempala_integrate (const int dimension, NT expected, NT exact, NT beta = 1.0, volumetype voltype = SOB) {

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef IsotropicQuadraticFunctor::parameters<NT> IsoParameters;

	IsoParameters params;
	params.alpha = (NT)2;
	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	HPOLYTOPE HP;
	HP = generate_cross <HPOLYTOPE> (dimension, false);
	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;
	NT integral_value = lovasz_vempala_integrate 
	  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, HP, x0, beta, voltype, 10, 0.1);
	
	test_values(integral_value, expected, exact);
}

template <typename NT>
void birkhoff_test_lovasz_vempala_integrate (const int dimension, NT expected, NT exact, NT beta = 1.0, volumetype voltype = SOB) {

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef IsotropicQuadraticFunctor::parameters<NT> IsoParameters;

	IsoParameters params;
	params.alpha = (NT)2;
	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	HPOLYTOPE HP;
	HP = generate_birkhoff <HPOLYTOPE> (dimension);
	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;
	NT integral_value = lovasz_vempala_integrate 
	  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, HP, x0, beta, voltype, 5, 0.1);
	
	test_values(integral_value, expected, exact);
}

template<typename NT>
void call_cubes_test_lovasz_vempala_integrate() {
	std::cout << "\nCUBE POLYTOPES\n";
	cubes_test_lovasz_vempala_integrate <NT> (1, 1.50, 1.49364, 1.0, CB);
	cubes_test_lovasz_vempala_integrate <NT> (2, 2.24, 2.23);
	cubes_test_lovasz_vempala_integrate <NT> (3, 3.34, 3.33225);
	cubes_test_lovasz_vempala_integrate <NT> (5, 7.50, 7.4325);
	cubes_test_lovasz_vempala_integrate <NT> (7, 16.59, 16.5852);
	// cubes_test_lovasz_vempala_integrate <NT> (10, 55.3, 55.26);
	// cubes_test_lovasz_vempala_integrate <NT> (15, 410.8554, 405);
	// cubes_test_lovasz_vempala_integrate <NT> (20, 3054.3493, 3050);
}

template<typename NT>
void call_simplices_test_lovasz_vempala_integrate() {
	std::cout << "\nSIMPLEX POLYTOPES\n";
	simplices_test_lovasz_vempala_integrate <NT> (1, 0.7392, 0.74711, 1.0, CB);
	simplices_test_lovasz_vempala_integrate <NT> (2, 0.365056, 0.366475);
	simplices_test_lovasz_vempala_integrate <NT> (3, 0.12797, 0.122755);
	simplices_test_lovasz_vempala_integrate <NT> (5, 0.00643434, 0.00659939);
	simplices_test_lovasz_vempala_integrate <NT> (7, 0.000165642, 0.000159737);
}

template<typename NT>
void call_prod_simplices_test_lovasz_vempala_integrate() {
	std::cout << "\nPRODUCT SIMPLEX POLYTOPES\n";
	prod_simplices_test_lovasz_vempala_integrate <NT> (1, 0.551034, 0.55361, 1.0, CB);
	prod_simplices_test_lovasz_vempala_integrate <NT> (2, 0.135389,  0.133205);
	prod_simplices_test_lovasz_vempala_integrate <NT> (3, 0.0155885, 0.0156483);
	prod_simplices_test_lovasz_vempala_integrate <NT> (5, 4.37432e-05, 4.36768e-05);
	prod_simplices_test_lovasz_vempala_integrate <NT> (7, 2.69889e-08, 2.70006e-08);
}

template<typename NT>
void call_cross_test_lovasz_vempala_integrate() {
	std::cout << "\nCROSS POLYTOPES\n";
	cross_test_lovasz_vempala_integrate <NT> (1, 1.50, 1.49364, 1.0, CB);
	cross_test_lovasz_vempala_integrate <NT> (2, 1.46701, 1.333333);
	cross_test_lovasz_vempala_integrate <NT> (3, 0.998024, 0.933333);
	cross_test_lovasz_vempala_integrate <NT> (5, 0.20977, 0.203174);
	cross_test_lovasz_vempala_integrate <NT> (7, 0.0208322, 0.020458);
}

template<typename NT>
void call_birkhoff_test_lovasz_vempala_integrate() {
	std::cout << "\nBIRKHOFF POLYTOPES\n";
	birkhoff_test_lovasz_vempala_integrate <NT> (2, 0.72793, 0.746982, 1.0, CB);
	birkhoff_test_lovasz_vempala_integrate <NT> (3, 0.0689502, 0.0682495);
	birkhoff_test_lovasz_vempala_integrate <NT> (4, 0.000425347, 0.000431226);
}

TEST_CASE("cubes") {
	call_cubes_test_lovasz_vempala_integrate<double>();
}

TEST_CASE("simplices") {
	call_simplices_test_lovasz_vempala_integrate<double>();
}

TEST_CASE("prod_simplices") {
	call_prod_simplices_test_lovasz_vempala_integrate<double>();
}

TEST_CASE("cross") {
	call_cross_test_lovasz_vempala_integrate<double>();
}

TEST_CASE("birkhoff") {
	call_birkhoff_test_lovasz_vempala_integrate<double>();
}


