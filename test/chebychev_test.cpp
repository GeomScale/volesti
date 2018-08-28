// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename NT, typename FilePath>
void cheb_test(FilePath f, NT expected, NT tolerance=0.0001)
{

    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open(f,std::ifstream::in);
    read_pointset(inp,Pin);
    int n = Pin[0][1]-1;
    HPolytope<Point> P;
    P.init(Pin);

    // Setup the parameters
    int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT,RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,true);

    //Compute chebychev ball//
    std::cout << "\n--- Testing Chebchev ball computation of " << f << std::endl;
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
    cheb_test<NT>("../data/cube10.ine", 1.0);
    cheb_test<NT>("../data/cube20.ine", 1.0);
    cheb_test<NT>("../data/cube30.ine", 1.0);
}

template <typename NT>
void call_test_cross() {
    cheb_test<NT>("../data/cross_10.ine", 0.316228);
}

template <typename NT>
void call_test_birk() {
    cheb_test<NT>("../data/birk3.ine", 0.207107);
    cheb_test<NT>("../data/birk4.ine", 0.122008);
    cheb_test<NT>("../data/birk5.ine", 0.0833333);
    cheb_test<NT>("../data/birk6.ine", 0.0618034);
}

template <typename NT>
void call_test_prod_simplex() {
    cheb_test<NT>("../data/prod_simplex_5_5.ine", 0.138197);
    cheb_test<NT>("../data/prod_simplex_10_10.ine", 0.0759747);
    cheb_test<NT>("../data/prod_simplex_15_15.ine", 0.0529858);
    cheb_test<NT>("../data/prod_simplex_20_20.ine", 0.0408628);
}

template <typename NT>
void call_test_simplex() {
    cheb_test<NT>("../data/simplex10.ine", 0.0759747);
    cheb_test<NT>("../data/simplex20.ine", 0.0408628);
    cheb_test<NT>("../data/simplex30.ine", 0.0281871);
    cheb_test<NT>("../data/simplex40.ine", 0.0215868);
    cheb_test<NT>("../data/simplex50.ine", 0.017522);
}

template <typename NT>
void call_test_skinny_cube() {
    cheb_test<NT>("../data/skinny_cube10.ine", 1.0);
    cheb_test<NT>("../data/skinny_cube20.ine", 1.0);
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
