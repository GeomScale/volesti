// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "polytope_generators.h"
#include <fstream>
#include <string>

template <class MT, class VT>
void create_txt(MT A, VT b, int kind, bool Vpoly) {

    int d = A.cols(), m = A.rows();
    std::string bar = "_";
    std::string ext = ".ext";
    std::string ine = ".ine";
    std::string name;

    std::ofstream outputFile;
    if(Vpoly) {
        if (kind == 0) {
            std::string poly = "zonotope";
            name = poly + bar + std::to_string(d) + bar + std::to_string(m) + ext;
            outputFile.open(name);
            outputFile<<name<<"\n";
            outputFile<<"Zonotpe\n";
        } else if (kind == 1) {
            std::string poly = "cube";
            name = poly + bar + std::to_string(d) + ext;
            outputFile.open(name);
            outputFile<<"cube_"<<d<<".ext\n";
            outputFile<<"V-representation\n";
        } else if (kind == 2) {
            std::string poly = "cross";
            name = poly + bar + std::to_string(d) + ext;
            outputFile.open(name);
            outputFile<<"cross_"<<d<<".ext\n";
            outputFile<<"V-representation\n";
        } else if (kind == 3) {
            std::string poly = "simplex";
            name = poly + bar + std::to_string(d) + ext;
            outputFile.open(name);
            outputFile<<"simplex_"<<d<<".ext\n";
            outputFile<<"V-representation\n";
        } else {
            return;
        }
    } else {
        if (kind == 1) {
            std::string poly = "cube";
            name = poly + bar + std::to_string(d) + ine;
            outputFile.open(name);
            outputFile<<"cube_"<<d<<".ine\n";
        } else if (kind == 2) {
            std::string poly = "cross";
            name = poly + bar + std::to_string(d) + ine;
            outputFile.open(name);
            outputFile<<"cross_"<<d<<".ine\n";
        } else if (kind == 3) {
            std::string poly = "simplex";
            name = poly + bar + std::to_string(d) + ine;
            outputFile.open(name);
            outputFile<<"simplex_"<<d<<".ine\n";
        } else if (kind == 4) {
            std::string poly = "prod_simplex";
            name = poly + bar + std::to_string(d/2) + bar + std::to_string(d/2) + ine;
            outputFile.open(name);
            outputFile<<"prod_simplex_"<<d/2<<"_"<<d/2<<".ine\n";
        } else if (kind == 5) {
            std::string poly = "skinny_cube";
            name = poly + bar + std::to_string(d) + ine;
            outputFile.open(name);
            outputFile<<"skinny_cube_"<<d<<".ine\n";
        } else {
            return;
        }
        outputFile<<"H-representation\n";
    }
    outputFile<<"begin\n";
    if(kind == 0){
        outputFile << " " << m << " " << d + 1 << " real\n";
    } else {
        outputFile << " " << m << " " << d + 1 << " integer\n";
    }

    for (int i = 0; i < m; ++i) {
        outputFile << " " << b(i);
        for (int j = 0; j < d; ++j) {
            outputFile << " " << A(i,j);
        }
        outputFile << "\n";
    }
    outputFile << "end\n";
    if (Vpoly) {
        outputFile << "hull\n";
        outputFile << "incidence";
    } else {
        outputFile << "input_incidence";
    }

    outputFile.close();
}

int main(const int argc, const char** argv) {
    //Deafault values
    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> Zonotope;

    Vpolytope VP;

    bool Hpoly = false, Vpoly = false, Zono = false,
            cube = false, cross = false, simplex = false,
                    prod_simplex = false, skinny_cube = false;
    int d = 0, m = 0, kind = -1;

    for(int i=1;i<argc;++i) {
        bool correct = false;

        if(!strcmp(argv[i],"-help")){
            std::cerr<<
                     "Usage:\n"<<
                     "-zonotope : generate a zonotope\n"<<
                     "-cube : generate a hypercube\n"<<
                     "-cross : generate a cross polytope\n"<<
                     "-simplex : generate a simplex\n"<<
                     "-prod_simplex : generate a product of two simplices\n"<<
                     "-skinny_cube : generate a skinny hypercube\n"<<
                     "-h : generate polytope in H-representation\n"<<
                     "-v : generate polytope in V-representation\n"<<
                     "-d : the dimension\n"<<
                     "-m : number of segments that generate the zonotope\n"<<
                     std::endl;
            return 0;
        }

        if(!strcmp(argv[i],"-zonotope")) {
            Zono = true;
            Vpoly = true;
            kind = 0;
            correct = true;
        }
        if(!strcmp(argv[i],"-cube")) {
            cube = true;
            kind = 1;
            correct = true;
        }
        if(!strcmp(argv[i],"-cross")) {
            cross = true;
            kind = 2;
            correct = true;
        }
        if(!strcmp(argv[i],"-simplex")) {
            simplex = true;
            kind = 3;
            correct = true;
        }
        if(!strcmp(argv[i],"-prod_simplex")) {
            prod_simplex = true;
            kind = 4;
            correct = true;
        }
        if(!strcmp(argv[i],"-skinny_cube")) {
            skinny_cube = true;
            kind = 5;
            correct = true;
        }

        if(!strcmp(argv[i],"-h")) {
            Hpoly = true;
            correct = true;
        }
        if(!strcmp(argv[i],"-v")) {
            Vpoly = true;
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

        if(correct==false) {
            std::cerr << "unknown parameters \'" << argv[i] <<
                      "\', try " << argv[0] << " --help" << std::endl;
            exit(-2);
        }

    }

    if (kind < 0 || d < 0) {
        std::cout<<"Wrong inputs, try -help"<<std::endl;
        exit(-1);
    }

    if (Zono) {
        if (m > 0) {
            Zonotope ZP = gen_zonotope<Zonotope, RNGType>(d, m);
            create_txt(ZP.get_mat(), ZP.get_vec(), kind, true);
        } else {
            std::cout << "Wrong inputs, try -help" << std::endl;
            exit(-1);
        }
    } else if (Hpoly) {
        Hpolytope HP;
        if (cube) {
            HP = gen_cube<Hpolytope>(d, false);
        } else if (cross) {
            HP = gen_cross<Hpolytope>(d, false);
        } else if (simplex) {
            HP = gen_simplex<Hpolytope>(d, false);
        } else if (prod_simplex) {
            HP = gen_prod_simplex<Hpolytope>(d);
        } else if (skinny_cube) {
            HP = gen_skinny_cube<Hpolytope>(d);
        } else {
            std::cout << "Wrong inputs, try -help" << std::endl;
            exit(-1);
        }
        create_txt(HP.get_mat(), HP.get_vec(), kind, false);
    } else if (Vpoly) {
        Vpolytope VP;
        if (cube) {
            VP = gen_cube<Vpolytope>(d, true);
        } else if (cross) {
            VP = gen_cross<Vpolytope>(d, true);
        } else if (simplex) {
            VP = gen_simplex<Vpolytope>(d, true);
        } else if (prod_simplex) {
            std::cout<<"No prod_simplex in V-representation can be generated, try -help"<<std::endl;
            exit(-1);
        } else if (skinny_cube) {
            std::cout<<"No skinny_cube in V-representation can be generated, try -help"<<std::endl;
            exit(-1);
        } else {
            std::cout<<"Wrong inputs, try -help"<<std::endl;
            exit(-1);
        }
        create_txt(VP.get_mat(), VP.get_vec(), kind, true);
    } else {
        std::cout<<"Wrong inputs, try -help"<<std::endl;
        exit(-1);
    }

    return 0;
}
