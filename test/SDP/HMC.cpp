//
// Created by panagiotis on 2/24/20.
//

#include "Eigen/Eigen"
#define VOLESTI_DEBUG
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
//#include "volume.h"
//#include "rotating.h"
#include "misc.h"
#include "linear_extensions.h"
//#include "cooling_balls.h"
//#include "cooling_hpoly.h"
//#include "sample_only.h"
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "matrices/EigenvaluesProblems.h"
#include "spectrahedron.h"
#include "LMI.h"
#include "matrices/EigenDenseMatrix.h"
#include "matrices/DenseProductMatrix.h"
#include "random_walks/HMC_RandomWalk.h"
#include "optimization/SDPA_FormatManager.h"
#include "random_walks/CoordinateDirectionsHitAndRun_RandomWalk.h"
#include "optimization/SimulatedAnnealing.h"


int main(int argc, char* argv[]) {


    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
    typedef HMC_RandomWalk<Point, MT, VT, RNGType > HMC_RandomWalk;
    typedef SimulatedAnnealing<Point, MT, VT> SA;

    std::ifstream inp;
    inp.open(argv[1],std::ifstream::in);
    LMI<NT, MT, VT> lmi;
    VT c;
    loadSDPAFormatFile(inp, lmi, c);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    SPECTRAHEDRON spectrahedron(lmi);
    lmi.print();
#ifdef SAMPLING
    unsigned int n = spectrahedron.dimension();

    RNGType rng(seed);

    int walkL=1, NN=30;
    NT T=0.1, diam=1;

//    SP.ComputeInnerBall(diam, radius);

    std::list<Point> randPoints;

    Point _c(c);
    HMC_RandomWalk::Settings settings(walkL, rng, _c, T, diam);
    HMC_RandomWalk hmc(settings);
    Point p(n);
    hmc.sample(spectrahedron, p, NN, randPoints);


    MT RetMat(n, NN);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = (*rpit).getCoefficients();

//    std::cout << RetMat;
    std::cout << "\n\n\n" << RetMat.transpose() * c << "\n";

#endif

    SA::Settings settings(0.001, 1,25,0.25);
    SA simulatedAnnealing(&spectrahedron, Point(c), settings);

    Point x;
    NT min = simulatedAnnealing.solve(x,true,-0.97226);

    std::cout << min << "\n";// << x.getCoefficients().transpose() << "\n";
    return 0;
}