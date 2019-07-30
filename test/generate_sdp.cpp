//
// Created by panagiotis on 28/5/2019.
//

#include "Eigen/Eigen"
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "convex_bodies/polytopes.h"
#include "polytope_generators.h"
#include <fstream>
#include <string>
#include "sdp_generator.h"
#include "sdp_problem.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef optimization::sdp_problem<Point> sdp_problem;




int main(const int argc, const char** argv) {
    //Deafault values
    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;

    bool file = false;
    std::ofstream os;

    int n = 0, m = 0;

    for(int i=1;i<argc;++i) {
        bool correct = false;

        if(!strcmp(argv[i],"-help")){
            std::cerr<<
                     "-n : the dimension of the objective function \n"<<
                     "-m : the dimension of the matrices \n"<<
                     "-f <filename> : the file to save the sdp" <<
                     std::endl;
            return 0;
        }



        if(!strcmp(argv[i],"-n")) {
            n = atof(argv[++i]);
            correct = true;
        }
        if(!strcmp(argv[i],"-m")) {
            m = atof(argv[++i]);
            correct = true;
        }

        if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--file")) {
            os.open(argv[++i]);
            file = true;
            correct = true;
        }

        if(correct==false) {
            std::cerr << "unknown parameters \'" << argv[i] <<
                      "\', try " << argv[0] << " --help" << std::endl;
            exit(-2);
        }

    }

    if (n <= 0 || !file ||  m<=0) {
        std::cout<<"Wrong inputs, try -help"<<std::endl;
        exit(-1);
    }

    sdp_problem sdp;
    sdp = generateSDP<Point>(n, m);

    sdp.writeSDPAFormatFile(os);
    os.close();



    return 0;
}
