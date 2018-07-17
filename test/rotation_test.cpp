#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"

long int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template <typename FilePath>
void rotation_test(FilePath f, double expected, double tolerance=0.1)//, double expected, double tolerance=0.1)
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
    std::cout << "\n--- Testing rotation of " << f << std::endl;
    double rot_val = rotating(P);
    //rot_val=std::abs(rot_val)
    std::pair<Point,double> CheBall = solveLP(P.get_matrix(), P.dimension());
    std::cout << "--- Testing volume of rotaded" << f << std::endl;
    double vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        vol += rot_val*volume(P,var,var,CheBall);
    }
    double error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    CHECK(error < tolerance);

}


TEST_CASE("rotated_skinny_cube") {
    rotation_test("../data/skinny_cube10.ine", 102400, 0.9);
    rotation_test("../data/skinny_cube20.ine", 104857600, 0.9);
}
