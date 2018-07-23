// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"

long int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename FilePath>
void test_volume(FilePath f, double expected, double tolerance=0.1)
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
    std::pair<Point,double> CheBall = P.chebyshev_center();

    // Estimate the volume
    std::cout << "--- Testing volume of " << f << std::endl;
    double vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        vol += volume(P,var,var,CheBall);
    }
    double error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    CHECK(error < tolerance);
}


TEST_CASE("cube") {
    test_volume("../data/cube10.ine", 1024);
    test_volume("../data/cube20.ine", 1048576);
    test_volume("../data/cube30.ine", 1073742000, 0.2);
}

TEST_CASE("cross") {
    test_volume("../data/cross_10.ine", 0.0002821869);
}

TEST_CASE("birkhoff") {
    test_volume("../data/birk3.ine", 0.125);
    test_volume("../data/birk4.ine", 0.000970018);
    test_volume("../data/birk5.ine", 0.000000225);
    test_volume("../data/birk6.ine", 0.0000000000009455459196, 0.5);
}

TEST_CASE("prod_simplex") {
    test_volume("../data/prod_simplex_5_5.ine", std::pow(1.0/factorial(5),2));
    test_volume("../data/prod_simplex_10_10.ine", std::pow(1.0/factorial(10),2));
    test_volume("../data/prod_simplex_15_15.ine", std::pow(1.0/factorial(15),2));
    test_volume("../data/prod_simplex_20_20.ine", std::pow(1.0/factorial(20),2));
}

TEST_CASE("simplex") {    
    test_volume("../data/simplex10.ine", 1.0/factorial(10));
    test_volume("../data/simplex20.ine", 1.0/factorial(20));
    test_volume("../data/simplex30.ine", 1.0/factorial(30));
    test_volume("../data/simplex40.ine", 1.0/factorial(40));
    test_volume("../data/simplex50.ine", 1.0/factorial(50));
}

TEST_CASE("skinny_cube") {    
    test_volume("../data/skinny_cube10.ine", 102400);
    test_volume("../data/skinny_cube20.ine", 104857600);
}
