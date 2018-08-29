// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"
#include <string>
#include <typeinfo>

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT, typename FilePath>
void test_CV_volume(FilePath f, NT expected, NT tolerance=0.2)
{
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open(f,std::ifstream::in);
    read_pointset(inp,Pin);
    int n = Pin[0][1]-1;
    VPolytope<Point, RNGType> VP;
    VP.init(Pin);

    // Setup the parameters
    int walk_len=1;
    int nexp=1, n_threads=1;
    NT e=0.2, err=0.0000000001;
    NT C=2.0,ratio,frac=0.1,delta=-1.0;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    int N = 500 * ((int) C) + ((int) (n * n / 2));
    int W = 4*n*n+500;
    ratio = 1.0-1.0/(NT(n));
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    //compute the chebychev ball
    std::pair<Point,NT> CheBall;


    // Estimate the volume
    std::cout << "--- Testing volume of " << f << std::endl;
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    NT vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        CheBall = VP.ComputeInnerBall();
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,rng,
                               urdist,urdist1,-1.0,false,false,false,false,false,false,true);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,false,
                                 false,false,false,false,false,false,true);
        vol += volume_gaussian_annealing(VP, var1, var2, CheBall);
        std::cout<<"vol = "<<vol<<std::endl;
    }
    NT error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
            CHECK(error < tolerance);
}

template <typename NT>
void call_test_cube(){
    test_CV_volume<NT>("../data/cube_3.ext", 8.0);
    test_CV_volume<NT>("../data/cube_4.ext", 16.0);
    test_CV_volume<NT>("../data/cube_5.ext", 32.0);
}

template <typename NT>
void call_test_cross(){
    test_CV_volume<NT>("../data/cross_4.ext", 0.6666667);
    if(typeid(NT)== typeid(double)) {
        test_CV_volume<NT>("../data/cross_5.ext", 0.26666667);
        test_CV_volume<NT>("../data/cross_6.ext", 0.08888889);
    }
}

template <typename NT>
void call_test_simplex() {
    test_CV_volume<NT>("../data/simplex5.ext", 1.0 / factorial(5.0));
    if(typeid(NT)== typeid(double)) {
        test_CV_volume<NT>("../data/simplex10.ext", 1.0 / factorial(10.0));
    }
}

TEST_CASE("cube") {
    call_test_cube<double>();
    call_test_cube<float>();
    call_test_cube<long double>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
    //call_test_cross<float>();
    call_test_cross<long double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
    //call_test_simplex<float>();
    call_test_simplex<long double>();
}
