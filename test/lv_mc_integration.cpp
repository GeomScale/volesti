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

template <typename NT>
void test_values (NT computed, NT expected) { //, NT exact) {
	std::cout << "----------------------------------------------------------Computed integration value = " << computed << std::endl;
	std::cout << "----------------------------------------------------------Expected integration value = " << expected << std::endl;
	// std::cout << "Exact integration value = " << exact << std::endl;
	// std::cout << "Relative error (expected) = " << std::abs((computed - expected)/expected) << std::endl;
	// std::cout << "Relative error (exact) = " << std::abs((computed - exact)/exact) << std::endl ;
	// CHECK(((std::abs((computed - expected)/expected) < 0.00001) || (std::abs((computed - exact)/exact) < 0.2)));
}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (1, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();;
	Point x0 = inner_ball.first;

	NT beta = 0.01;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, CB, 5, 0.1);
	
	test_values(integral_value, 1.49364);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate2() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (2, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();;
	Point x0 = inner_ball.first;

	NT beta = 0.01;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 2.23);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate3() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (3, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();;
	Point x0 = inner_ball.first;

	NT beta = 0.01;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 3.33);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate4() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (4, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 0.01;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 4.977);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate5() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (5, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 1.0;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 7.434);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate6() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (6, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 1.0;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 11.1039);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate7() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (7, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 1.0;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 16.58);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate8() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (8, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 1.0;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 24.77);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate9() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (9, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 1.0;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 37.001);

}

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate10() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (10, false);

	IsotropicQuadraticFunctor::parameters<NT> params;
	params.alpha = (NT) 2;

	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	unsigned int n = HP.dimension();
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

	NT beta = 1.0;

	NT integral_value = lovasz_vempala_integrate 
	  <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, BilliardWalk, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, beta, SOB, 5, 0.1);
	
	test_values(integral_value, 55.266);

}


TEST_CASE("iso") {
	call_cubes_test_lovasz_vempala_integrate<double>();
}

TEST_CASE("iso2") {
	call_cubes_test_lovasz_vempala_integrate2<double>();
}

TEST_CASE("iso3") {
	call_cubes_test_lovasz_vempala_integrate3<double>();
}

TEST_CASE("iso4") {
	call_cubes_test_lovasz_vempala_integrate4<double>();
}

TEST_CASE("iso5") {
	call_cubes_test_lovasz_vempala_integrate5<double>();
}

// TEST_CASE("iso6") {
// 	call_cubes_test_lovasz_vempala_integrate6<double>();
// }

// TEST_CASE("iso7") {
// 	call_cubes_test_lovasz_vempala_integrate7<double>();
// }

// TEST_CASE("iso8") {
// 	call_cubes_test_lovasz_vempala_integrate8<double>();
// }

// TEST_CASE("iso9") {
// 	call_cubes_test_lovasz_vempala_integrate9<double>();
// }

// TEST_CASE("iso10") {
// 	call_cubes_test_lovasz_vempala_integrate10<double>();
// }