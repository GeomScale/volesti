// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#define VOLESTI_DEBUG
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "rotating.h"
#include "misc.h"
#include "linear_extensions.h"
#include "cooling_balls.h"
#include "cooling_hpoly.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "Eigen/Eigen"
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "rotating.h"
#include "misc.h"
#include "doctest.h"
#include <unistd.h>
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "cooling_hpoly.h"
#include "cooling_balls.h"
#include "misc.h"
#include "polytope_generators.h"
#include <typeinfo>
#include "sample_only.h"
#include "exact_vols.h"

template <typename NT>
NT factorial(NT n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT, class RNGType, class Polytope>
void test_volume(Polytope &HP, NT expected, NT tolerance=0.1)
{

    typedef typename Polytope::PolytopePoint Point;

    // Setup the parameters
    int n = HP.dimension();
    int walk_len=1;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::pair<Point, NT> InnerBall;
    InnerBall = HP.ComputeInnerBall();
    NT diameter;
    HP.comp_diam(diameter, InnerBall.second);

    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,InnerBall.second,diameter,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,true);

    //Compute chebychev ball//
    std::pair<Point,NT> CheBall;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    NT vol = 0;
    NT ball_radius=0.0, lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2, round_val = 1.0;
    NT C=2.0,ratio,frac=0.1,delta=-1.0;
    vars_ban <NT> var_ban(lb, ub, p, rmax, alpha, 150, 125, 10, false);

    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        CheBall = HP.ComputeInnerBall();
        vol += vol_cooling_balls(HP, var, var_ban, InnerBall);
    }
    NT error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    CHECK(error < tolerance);
}


template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    test_volume<NT, RNGType>(P, 1024.0);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    test_volume<NT, RNGType>(P, 1048576.0);

    std::cout << "--- Testing volume of H-cube30" << std::endl;
    P = gen_cube<Hpolytope>(30, false);
    test_volume<NT, RNGType>(P, 1073742000.0, 0.2);
}

template <typename NT>
void call_test_cross(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;

    std::cout << "--- Testing volume of H-cross10" << std::endl;
    Hpolytope P = gen_cross<Hpolytope>(10, false);
    test_volume<NT, RNGType>(P, 0.0002821869);
}

template <typename NT>
void call_test_birk() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-birk3" << std::endl;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open("../R-proj/inst/extdata/birk3.ine",std::ifstream::in);
    read_pointset(inp,Pin);
    P.init(Pin);
    test_volume<NT, RNGType>(P, 0.125);

    std::cout << "--- Testing volume of H-birk4" << std::endl;
    std::ifstream inp2;
    std::vector<std::vector<NT> > Pin2;
    inp2.open("../R-proj/inst/extdata/birk4.ine",std::ifstream::in);
    read_pointset(inp2,Pin2);
    P.init(Pin2);
    test_volume<NT, RNGType>(P, 0.000970018);

    std::cout << "--- Testing volume of H-birk5" << std::endl;
    std::ifstream inp3;
    std::vector<std::vector<NT> > Pin3;
    inp3.open("../R-proj/inst/extdata/birk5.ine",std::ifstream::in);
    read_pointset(inp3,Pin3);
    P.init(Pin3);
    test_volume<NT, RNGType>(P, 0.000000225);

    std::cout << "--- Testing volume of H-birk6" << std::endl;
    std::ifstream inp4;
    std::vector<std::vector<NT> > Pin4;
    inp4.open("../R-proj/inst/extdata/birk6.ine",std::ifstream::in);
    read_pointset(inp4,Pin4);
    P.init(Pin4);
    test_volume<NT, RNGType>(P, 0.0000000000009455459196, 0.5);
}

template <typename NT>
void call_test_prod_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-prod_simplex5" << std::endl;
    P = gen_prod_simplex<Hpolytope>(5);
    test_volume<NT, RNGType>(P, std::pow(1.0 / factorial(5.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
    P = gen_prod_simplex<Hpolytope>(10);
    test_volume<NT, RNGType>(P, std::pow(1.0 / factorial(10.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
    P = gen_prod_simplex<Hpolytope>(15);
    test_volume<NT, RNGType>(P, std::pow(1.0 / factorial(15.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex20" << std::endl;
    P = gen_prod_simplex<Hpolytope>(20);
    test_volume<NT, RNGType>(P, std::pow(1.0 / factorial(20.0), 2));
}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-simplex10" << std::endl;
    P = gen_simplex<Hpolytope>(10, false);
    test_volume<NT, RNGType>(P, 1.0 / factorial(10.0));

    std::cout << "--- Testing volume of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    test_volume<NT, RNGType>(P, 1.0 / factorial(20.0));

    std::cout << "--- Testing volume of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    test_volume<NT, RNGType>(P, 1.0 / factorial(30.0));

    std::cout << "--- Testing volume of H-simplex40" << std::endl;
    P = gen_simplex<Hpolytope>(40, false);
    test_volume<NT, RNGType>(P, 1.0 / factorial(40.0));

    //std::cout << "--- Testing volume of H-simplex50" << std::endl;
    //P = gen_simplex<Hpolytope>(50, false);
    //test_volume<NT, RNGType>(P, 1.0 / factorial(50.0));
}

template <typename NT>
void call_test_skinny_cube() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-skinny_cube10" << std::endl;
    P = gen_skinny_cube<Hpolytope>(10);
    test_volume<NT, RNGType>(P, 102400.0);

    //std::cout << "--- Testing volume of H-skinny_cube20" << std::endl;
    //P = gen_skinny_cube<Hpolytope>(20);
    //test_volume<NT, RNGType>(P, 104857600.0);
}


TEST_CASE("cube") {
    call_test_cube<double>();
    //call_test_cube<float>();
    //call_test_cube<long double>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
    //call_test_cross<float>();
    //call_test_cross<long double>();
}

TEST_CASE("birk") {
    call_test_birk<double>();
    //call_test_birk<float>();
    //call_test_birk<long double>();
}

TEST_CASE("prod_simplex") {
    call_test_prod_simplex<double>();
    //call_test_prod_simplex<float>();
    //call_test_prod_simplex<long double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
    //call_test_simplex<float>();
    //call_test_simplex<long double>();
}

TEST_CASE("skinny_cube") {
    call_test_skinny_cube<double>();
    //call_test_skinny_cube<float>();
    //call_test_skinny_cube<long double>();
}
