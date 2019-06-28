
#include "Eigen/Eigen"

#define VOLESTI_DEBUG


#include <fstream>
#include "volume.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "solve_lp.h"

#include "cutting_plane_sdp.h"
#include "sdp_problem.h"
#include <chrono>
#include "Eigen"
#include <iomanip>

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef boost::mt19937 RNGType;
typedef HPolytope<Point> Hpolytope;
//typedef optimization::lp_problem<Point, NT> lp_problem;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

extern int STEPS; //TODO delete this

void printHelpMessage();

int main(const int argc, const char **argv) {

    return 0;
}


