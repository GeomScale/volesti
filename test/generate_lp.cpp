// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Panagiotis Repouskos, as part of Google Summer of Code 2019 program.

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
#include "lp_problem.h"
#include "lp_generator.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef optimization::lp_problem<Point, NT> lp_problem;

typedef enum kind {
    cube, cross, simplex, prod_simplex, skinny_cube, zonotope, hpoly
} Kind;



int main(const int argc, const char** argv) {
    //Deafault values
    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> Zonotope;

    bool file = false;
    std::ofstream os;
    Vpolytope VP;

    int d = 0, m = 0;
    Kind kind = cube;

    for(int i=1;i<argc;++i) {
        bool correct = false;

        if(!strcmp(argv[i],"-help")){
            std::cerr<<
//                     "Usage:\n"<<
//                     "-zonotope : generate a zonotope\n"<<
//                     "-cube : generate a hypercube\n"<<
//                     "-cross : generate a cross polytope\n"<<
//                     "-simplex : generate a simplex\n"<<
//                     "-prod_simplex : generate a product of two simplices\n"<<
//                     "-skinny_cube : generate a skinny hypercube\n"<<
//                     "-h : generate polytope in H-representation\n"<<
//                     "-v : generate polytope in V-representation\n"<<
                     "-d : the dimension\n"<<
//                     "-m : number of segments that generate the zonotope\n"<<
                      "-f <filename> : the file to save the lp" <<
                     std::endl;
            return 0;
        }

        if(!strcmp(argv[i],"-zonotope")) {
            kind = zonotope;
            correct = true;
        }
        if(!strcmp(argv[i],"-cube")) {
            kind = cube;
            correct = true;
        }
        if(!strcmp(argv[i],"-cross")) {
            kind = cross;
            correct = true;
        }
        if(!strcmp(argv[i],"-simplex")) {
            kind = simplex;
            correct = true;
        }
        if(!strcmp(argv[i],"-prod_simplex")) {
            kind = prod_simplex;
            correct = true;
        }
        if(!strcmp(argv[i],"-skinny_cube")) {
            kind = skinny_cube;
            correct = true;
        }

        if(!strcmp(argv[i],"-h")) {
            kind = hpoly;
            correct = true;
        }

        if(!strcmp(argv[i],"-d")) {
            d = atof(argv[++i]);
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

    if (d <= 0 || !file || (kind == zonotope && m<=0) || (kind == hpoly && m<=0)) {
        std::cout<<"Wrong inputs, try -help"<<std::endl;
        exit(-1);
    }

    lp_problem lp;

    switch (kind) {
        case cube:
            lp= generate_lp_cube<Point, NT>(d);
            break;
        case cross:
            lp= generate_lp_cross<Point, NT>(d);
            break;
        case simplex:
            lp= generate_lp_simplex<Point, NT>(d);
            break;
        case prod_simplex:
            lp= generate_lp_prod_simplex<Point, NT>(d);
            break;
        case skinny_cube:
            lp= generate_lp_skinny_cube<Point, NT>(d);
            break;
        case zonotope:
            lp= generate_lp_zonotope<Point, NT, RNGType>(d, m);
            break;
        case hpoly:
            lp= generate_lp<Point, NT, RNGType>(d, m);
            break;
    }

    lp.saveToFile(os);
    os.close();


    return 0;
}
