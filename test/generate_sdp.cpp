//
// Created by panagiotis on 28/5/2019.
//

#include "Eigen/Eigen"
#include <chrono>
#include <fstream>
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include "boost/random.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/uniform_real_distribution.hpp"
#include "vars.h"
#include "samplers.h"
#include "rounding.h"
#include "sample_only.h"
#include "sdp_generator.h"
#include "spectrahedron.h"


int main(const int argc, const char** argv) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT, NT, VT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <MT, VT> lmi;
    typedef Spectrahedron <lmi, Point> spectaedro;

    bool file = false;
    //std::ofstream os;

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
            //os.open(argv[++i]);
            //file = true;
            correct = true;
        }

        if(correct==false) {
            std::cerr << "unknown parameters \'" << argv[i] <<
                      "\', try " << argv[0] << " --help" << std::endl;
            exit(-2);
        }

    }

    if (n <= 0 ||  m<=0) {
        std::cout<<"Wrong inputs, try -help"<<std::endl;
        exit(-1);
    }

    spectaedro SP;//, SP2;
    SP = generateSDP2<lmi, spectaedro, Point>(n, m);

    Point c = get_direction<RNGType, Point, NT>(n);

    std::filebuf fb;
    std::string bar = "_";
    std::string txt = ".txt";
    fb.open("sdp_prob"+bar+std::to_string(n)+bar+std::to_string(m)+txt, std::ios::out);
    std::ostream os(&fb);
    writeSDPAFormatFile<MT>(os, SP.getLMI(), c.get_coefficients());

    //sdp_problem sdp;
   // sdp = generateSDP<Point>(n, m);

    //sdp.writeSDPAFormatFile(os);
    //os.close();



    return 0;
}
