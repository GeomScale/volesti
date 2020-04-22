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

//    if (Zono) {
//        if (m > 0) {
//            Zonotope ZP = gen_zonotope<Zonotope, RNGType>(d, m);
//            create_txt(ZP.get_mat(), ZP.get_vec(), kind, true);
//        } else {
//            std::cout << "Wrong inputs, try -help" << std::endl;
//            exit(-1);
//        }
//    } else if (Hpoly) {
//        Hpolytope HP;
//        if (cube) {
//            HP = gen_cube<Hpolytope>(d, false);
//        } else if (cross) {
//            HP = gen_cross<Hpolytope>(d, false);
//        } else if (simplex) {
//            HP = gen_simplex<Hpolytope>(d, false);
//        } else if (prod_simplex) {
//            HP = gen_prod_simplex<Hpolytope>(d);
//        } else if (skinny_cube) {
//            HP = gen_skinny_cube<Hpolytope>(d);
//        } else {
//            std::cout << "Wrong inputs, try -help" << std::endl;
//            exit(-1);
//        }
//        create_txt(HP.get_mat(), HP.get_vec(), kind, false);
//    } else if (Vpoly) {
//        Vpolytope VP;
//        if (cube) {
//            VP = gen_cube<Vpolytope>(d, true);
//        } else if (cross) {
//            VP = gen_cross<Vpolytope>(d, true);
//        } else if (simplex) {
//            VP = gen_simplex<Vpolytope>(d, true);
//        } else if (prod_simplex) {
//            std::cout<<"No prod_simplex in V-representation can be generated, try -help"<<std::endl;
//            exit(-1);
//        } else if (skinny_cube) {
//            std::cout<<"No skinny_cube in V-representation can be generated, try -help"<<std::endl;
//            exit(-1);
//        } else {
//            std::cout<<"Wrong inputs, try -help"<<std::endl;
//            exit(-1);
//        }
//        create_txt(VP.get_mat(), VP.get_vec(), kind, true);
//    } else {
//        std::cout<<"Wrong inputs, try -help"<<std::endl;
//        exit(-1);
//    }

    return 0;
}
