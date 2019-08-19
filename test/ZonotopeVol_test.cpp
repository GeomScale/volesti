// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"
#include "polytope_generators.h"
#include "exact_vols.h"
#include <typeinfo>

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT>
void test_zono_volume(int n, int m, NT tolerance = 0.15)
{
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef Zonotope<Point> Zonotope;
    Zonotope ZP = gen_zonotope<Zonotope, RNGType>(n, m);

    // Setup the parameters
    int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001, round_value = 1.0;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    std::pair<NT,NT> res_round;

    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                          urdist,urdist1,-1.0,false,false,false,false,false,false,true,false,false);

    //Compute chebychev ball//
    std::pair<Point,NT> CheBall;
    CheBall = ZP.ComputeInnerBall();

    // Estimate the volume
    std::cout << "--- Testing volume of Zonotope in dimension: " << n <<" and number of generators: "<< m << std::endl;
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    NT vol_exact = exact_zonotope_vol<NT>(ZP);
    res_round = rounding_min_ellipsoid(ZP, CheBall, var);
    round_value = round_value * res_round.first;
    NT vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        CheBall = ZP.ComputeInnerBall();
        vol += round_value * volume(ZP,var,CheBall);
    }
    NT error = std::abs(((vol/num_of_exp)-vol_exact))/vol_exact;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << vol_exact << std::endl;
            CHECK(error < tolerance);
}


template <typename NT>
void call_test(int n, int m){
    test_zono_volume<NT>(n, m);
}


TEST_CASE("4_dimensional") {
    call_test<double>(4, 8);
    //call_test<float>(4,8);
    //call_test<long double>(4,8);
}
