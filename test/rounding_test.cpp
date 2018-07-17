#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"

long int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename FilePath>
void rounding_test(FilePath f)//, double expected, double tolerance=0.1)
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

    vars var(rnum,n,walk_len,n_threads,err,e,0,0,0,rng,
             urdist,urdist1,false,false,false,false,false,true);

    //Compute chebychev ball//
    std::cout << "\n--- Testing rounding of " << f << std::endl;
    std::pair<Point,double> CheBall = solveLP(P.get_matrix(), P.dimension());
    Point c=CheBall.first;
    NT radius=CheBall.second;
    NT round_value;
    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    round_value = rounding_min_ellipsoid(P,c,radius,var);
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;

    std::cout<<"\nround value is: "<<round_value<<std::endl;
    std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;

    CHECK(!std::isnan(round_value));
    CHECK(!std::isinf(round_value));
}

TEST_CASE("round_cube") {
    rounding_test("../data/cube10.ine");
    rounding_test("../data/cube20.ine");
    rounding_test("../data/cube30.ine");
}

TEST_CASE("round_cross") {
    rounding_test("../data/cross_10.ine");
}

TEST_CASE("round_birkhoff") {
    rounding_test("../data/birk3.ine");
    rounding_test("../data/birk4.ine");
    rounding_test("../data/birk5.ine");
    rounding_test("../data/birk6.ine");
}

TEST_CASE("round_prod_simplex") {
    rounding_test("../data/prod_simplex_5_5.ine");
    rounding_test("../data/prod_simplex_10_10.ine");
    rounding_test("../data/prod_simplex_15_15.ine");
    rounding_test("../data/prod_simplex_20_20.ine");
}

TEST_CASE("round_simplex") {
    rounding_test("../data/simplex10.ine");
    rounding_test("../data/simplex20.ine");
    rounding_test("../data/simplex30.ine");
    rounding_test("../data/simplex40.ine");
    rounding_test("../data/simplex50.ine");
}

TEST_CASE("round_skinny_cube") {
    rounding_test("../data/skinny_cube10.ine");
    rounding_test("../data/skinny_cube20.ine");
}