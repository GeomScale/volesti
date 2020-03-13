// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "misc.h"
#include "known_polytope_generators.h"

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename NT, class RNGType, class Polytope>
void cheb_test(Polytope &P, NT expected, NT tolerance=0.0001)
{

    typedef typename Polytope::PolytopePoint Point;

    // Setup the parameters
    int n = P.dimension();
    int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT,RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,true, false,false);

    //Compute chebychev ball//
    //std::cout << "\n--- Testing Chebchev ball computation of " << f << std::endl;
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    std::pair<Point,NT> CheBall = P.ComputeInnerBall();
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;


    std::cout<<"\nradius is: "<<CheBall.second<<std::endl;
    std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;

    NT error = std::abs(CheBall.second-expected)/expected;

    CHECK(error < tolerance);

}

template <typename NT>
void call_test_cube() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    cheb_test<NT, RNGType>(P, 1.0);

    std::cout << "\n--- Testing Chebchev ball computation of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    cheb_test<NT, RNGType>(P, 1.0);

    std::cout << "\n--- Testing Chebchev ball computation of H-cube30" << std::endl;
    P = gen_cube<Hpolytope>(30, false);
    cheb_test<NT, RNGType>(P, 1.0);

}

template <typename NT>
void call_test_cross() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;

    std::cout << "\n--- Testing Chebchev ball computation of H-cross10" << std::endl;
    Hpolytope P = gen_cross<Hpolytope>(10, false);
    cheb_test<NT, RNGType>(P, 0.316228);

}

template <typename NT>
void call_test_birk() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-birk3" << std::endl;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open("../R-proj/inst/extdata/birk3.ine",std::ifstream::in);
    read_pointset(inp,Pin);
    P.init(Pin);
    cheb_test<NT, RNGType>(P, 0.207107);

    std::cout << "\n--- Testing Chebchev ball computation of H-birk4" << std::endl;
    std::ifstream inp2;
    std::vector<std::vector<NT> > Pin2;
    inp2.open("../R-proj/inst/extdata/birk4.ine",std::ifstream::in);
    read_pointset(inp2,Pin2);
    P.init(Pin2);
    cheb_test<NT, RNGType>(P, 0.122008);

    std::cout << "\n--- Testing Chebchev ball computation of H-birk5" << std::endl;
    std::ifstream inp3;
    std::vector<std::vector<NT> > Pin3;
    inp3.open("../R-proj/inst/extdata/birk5.ine",std::ifstream::in);
    read_pointset(inp3,Pin3);
    P.init(Pin3);
    cheb_test<NT, RNGType>(P, 0.0833333);

    std::cout << "\n--- Testing Chebchev ball computation of H-birk6" << std::endl;
    std::ifstream inp4;
    std::vector<std::vector<NT> > Pin4;
    inp4.open("../R-proj/inst/extdata/birk6.ine",std::ifstream::in);
    read_pointset(inp4,Pin4);
    P.init(Pin4);
    cheb_test<NT, RNGType>(P, 0.0618034);

}

template <typename NT>
void call_test_prod_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-prod_simplex5" << std::endl;
    P = gen_prod_simplex<Hpolytope>(5);
    cheb_test<NT, RNGType>(P, 0.138197);

    std::cout << "\n--- Testing Chebchev ball computation of H-prod_simplex10" << std::endl;
    P = gen_prod_simplex<Hpolytope>(10);
    cheb_test<NT, RNGType>(P, 0.0759747);

    std::cout << "\n--- Testing Chebchev ball computation of H-prod_simplex15" << std::endl;
    P = gen_prod_simplex<Hpolytope>(15);
    cheb_test<NT, RNGType>(P, 0.0529858);

    std::cout << "\n--- Testing Chebchev ball computation of H-prod_simplex20" << std::endl;
    P = gen_prod_simplex<Hpolytope>(20);
    cheb_test<NT, RNGType>(P, 0.0408628);

}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex10" << std::endl;
    P = gen_simplex<Hpolytope>(10, false);
    cheb_test<NT, RNGType>(P, 0.0759747);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    cheb_test<NT, RNGType>(P, 0.0408628);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    cheb_test<NT, RNGType>(P, 0.0281871);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex40" << std::endl;
    P = gen_simplex<Hpolytope>(40, false);
    cheb_test<NT, RNGType>(P, 0.0215868);

    std::cout << "\n--- Testing Chebchev ball computation of H-simplex50" << std::endl;
    P = gen_simplex<Hpolytope>(50, false);
    cheb_test<NT, RNGType>(P, 0.017522);

}

template <typename NT>
void call_test_skinny_cube() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing Chebchev ball computation of H-skinny_cube10" << std::endl;
    P = gen_skinny_cube<Hpolytope>(10);
    cheb_test<NT, RNGType>(P, 1.0);

    std::cout << "\n--- Testing Chebchev ball computation of H-skinny_cube20" << std::endl;
    P = gen_skinny_cube<Hpolytope>(20);
    cheb_test<NT, RNGType>(P, 1.0);

}



TEST_CASE("cheb_cube") {
    call_test_cube<double>();
    call_test_cube<float>();
    call_test_cube<long double>();
}

TEST_CASE("cheb_cross") {
    call_test_cross<double>();
    call_test_cross<float>();
    call_test_cross<long double>();
}

TEST_CASE("cheb_birk") {
    call_test_birk<double>();
    call_test_birk<float>();
    call_test_birk<long double>();
}

TEST_CASE("cheb_prod_simplex") {
    call_test_prod_simplex<double>();
    call_test_prod_simplex<float>();
    call_test_prod_simplex<long double>();
}

TEST_CASE("cheb_simplex") {
    call_test_simplex<double>();
    call_test_simplex<float>();
    call_test_simplex<long double>();
}

TEST_CASE("cheb_skinny_cube") {
    call_test_skinny_cube<double>();
    call_test_skinny_cube<float>();
    call_test_skinny_cube<long double>();
}
