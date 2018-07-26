// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"

long int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename FilePath>
void rounding_test(FilePath f, bool rot, double expected, double tolerance=0.1)//, double expected, double tolerance=0.1)
{
    std::ifstream inp;
    std::vector<std::vector<double> > Pin;
    inp.open(f,std::ifstream::in);
    read_pointset(inp,Pin);
    int n = Pin[0][1]-1;
    Polytope<double> P;
    P.init(Pin);

    // Setup the parameters
    double walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    double e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,true);

    //Compute chebychev ball//
    std::cout << "\n--- Testing rounding of " << f << std::endl;
    double rot_val;
    if(rot){
        std::cout << "\n--- Rotation is ON "<< std::endl;
        rot_val = rotating(P);
        std::cout << "Rotation value = "<<rot_val<<std::endl;
    }
    std::pair<Point,double> CheBall;// = solveLP(P.get_matrix(), P.dimension());
    Point c;//=CheBall.first;
    NT radius;//=CheBall.second;
    NT round_value=1.0, ratio1,ratio2;
    std::pair<NT,NT> res_round;
    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    int count=1;
    CheBall = P.chebyshev_center();
    //c=CheBall.first;
    //radius=CheBall.second;
    res_round = rounding_min_ellipsoid(P, CheBall, var);
    round_value = round_value * res_round.first;
    ratio2 = res_round.second;
    ratio1 = 0.0;
    //std::cout<<ratio1<<" "<<ratio2<<std::endl;
    while(ratio2>ratio1 && count<=4) {
        CheBall = P.chebyshev_center();
        //c=CheBall.first;
        //radius=CheBall.second;
        res_round = rounding_min_ellipsoid(P, CheBall, var);
        round_value = round_value * res_round.first;
        ratio1=ratio2;
        ratio2 = res_round.second;
        //std::cout<<ratio1<<" "<<ratio2<<std::endl;
        count++;
    }
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout<<"\nround value is: "<<round_value<<std::endl;
    std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
    CheBall = P.chebyshev_center();
    //c=CheBall.first;
    //radius=CheBall.second;

    double vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        vol += round_value*volume(P,var,var,CheBall);
    }
    double error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    CHECK(error < tolerance);



}

TEST_CASE("round_rot_skinny_cube") {
    rounding_test("../data/skinny_cube10.ine", true, 102400);
    rounding_test("../data/skinny_cube20.ine", true, 104857600, 0.2);
}

TEST_CASE("round_skinny_cube") {
    rounding_test("../data/skinny_cube10.ine", false, 102400);
    rounding_test("../data/skinny_cube20.ine", false, 104857600);
}
