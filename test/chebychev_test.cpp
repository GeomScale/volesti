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
    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    std::pair<Point,NT> CheBall = P.ComputeInnerBall();
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;


    std::cout<<"\nradius is: "<<CheBall.second<<std::endl;
    std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;

    NT error = std::abs(CheBall.second-expected)/expected;

    CHECK(error < tolerance);

}

template <typename NT, class string>
void call_tests(string test) {

    string strcube("cheb_cube");
    string strcross("cheb_cross");
    string strbirk("cheb_birk");
    string strprod_simplex("cheb_prod_simplex");
    string strsimplex("cheb_simplex");
    string strskinny_cube("skinny_cube");

    if(test.compare(strcube) == 0) {
        cheb_test<NT>("../data/cube10.ine", 1.0);
        cheb_test<NT>("../data/cube20.ine", 1.0);
        cheb_test<NT>("../data/cube30.ine", 1.0);
    }

    if(test.compare(strcross) == 0) {
        cheb_test<NT>("../data/cross_10.ine", 0.316228);
    }

    if(test.compare(strbirk) == 0) {
        cheb_test<NT>("../data/birk3.ine", 0.207107);
        cheb_test<NT>("../data/birk4.ine", 0.122008);
        cheb_test<NT>("../data/birk5.ine", 0.0833333);
        cheb_test<NT>("../data/birk6.ine", 0.0618034);
    }

    if(test.compare(strprod_simplex) == 0) {
        cheb_test<NT>("../data/prod_simplex_5_5.ine", 0.138197);
        cheb_test<NT>("../data/prod_simplex_10_10.ine", 0.0759747);
        cheb_test<NT>("../data/prod_simplex_15_15.ine", 0.0529858);
        cheb_test<NT>("../data/prod_simplex_20_20.ine", 0.0408628);
    }

    if(test.compare(strsimplex) == 0) {
        cheb_test<NT>("../data/simplex10.ine", 0.0759747);
        cheb_test<NT>("../data/simplex20.ine", 0.0408628);
        cheb_test<NT>("../data/simplex30.ine", 0.0281871);
        cheb_test<NT>("../data/simplex40.ine", 0.0215868);
        cheb_test<NT>("../data/simplex50.ine", 0.017522);
    }

    if(test.compare(strskinny_cube) == 0) {
        cheb_test<NT>("../data/skinny_cube10.ine", 1.0);
        cheb_test<NT>("../data/skinny_cube20.ine", 1.0);
    }
}



TEST_CASE("cheb_cube") {
    std::string strtest("cheb_cube");
    call_tests<double>(strtest);
    call_tests<float>(strtest);
    call_tests<long double>(strtest);
}

TEST_CASE("cheb_cross") {
    std::string strtest("cheb_cross");
    call_tests<double>(strtest);
    call_tests<float>(strtest);
    call_tests<long double>(strtest);
}

TEST_CASE("cheb_birk") {
    std::string strtest("cheb_birk");
    call_tests<double>(strtest);
    call_tests<float>(strtest);
    call_tests<long double>(strtest);
}

TEST_CASE("cheb_prod_simplex") {
    std::string strtest("cheb_prod_simplex");
    call_tests<double>(strtest);
    call_tests<float>(strtest);
    call_tests<long double>(strtest);
}

TEST_CASE("cheb_simplex") {
    std::string strtest("cheb_simplex");
    call_tests<double>(strtest);
    call_tests<float>(strtest);
    call_tests<long double>(strtest);
}

TEST_CASE("cheb_skinny_cube") {
    std::string strtest("cheb_skinny_cube");
    call_tests<double>(strtest);
    call_tests<float>(strtest);
    call_tests<long double>(strtest);
}
