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
#include "vpolyoracles.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "known_polytope_generators.h"
#include "h_polytopes_gen.h"
#include "v_polytopes_gen.h"
#include "z_polytopes_gen.h"
#include <fstream>
#include <string>

template <class MT, class VT>
void create_txt(MT A, VT b, int kind, bool Vpoly, bool Zono) {

    int d = A.cols(), m = A.rows();
    std::string bar = "_", ext = ".ext", ine = ".ine", name;

    std::ofstream outputFile;

    if (Zono) {
        name = "zonotope" + bar + std::to_string(d) + bar + std::to_string(m) + ext;
        outputFile.open(name);
        outputFile << name << "\n";
        outputFile << "Zonotope\n";
    } else if (Vpoly) {
        switch (kind) {
            case 1:
                name = "cube" + bar + std::to_string(d) + ext;
                break;
            case 2:
                name = "cross" + bar + std::to_string(d) + ext;
                break;
            case 3:
                name = "simplex" + bar + std::to_string(d) + ext;
                break;
            case 4:
                name = "rvc" + bar + std::to_string(d) + bar + std::to_string(m) + ext;
                break;
            case 5:
                name = "rvs" + bar + std::to_string(d) + bar + std::to_string(m) + ext;
                break;
            default:
                return;
        }
        outputFile.open(name);
        outputFile << name << "\n";
        outputFile << "V-representation\n";
    } else {
        switch (kind) {
            case 1:
                name = "cube" + bar + std::to_string(d) + ine;
                break;
            case 2:
                name = "cross" + bar + std::to_string(d) + ine;
                break;
            case 3:
                name = "simplex" + bar + std::to_string(d) + ine;
                break;
            case 4:
                name = "prod_simplex" + bar + std::to_string(d / 2) + bar + std::to_string(d / 2) + ine;
                break;
            case 5:
                name = "skinny_cube" + bar + std::to_string(d) + ine;
                break;
            case 6:
                name = "rhs" + bar + std::to_string(d) + bar + std::to_string(m) + ine;
                break;
            default:
                return;
        }
        outputFile.open(name);
        outputFile << name << "\n";
        outputFile << "H-representation\n";
    }
    outputFile << "begin\n";
    if (kind == 0) {
        outputFile << " " << m << " " << d + 1 << " real\n";
    } else {
        outputFile << " " << m << " " << d + 1 << " integer\n";
    }

    for (int i = 0; i < m; ++i) {
        outputFile << " " << b(i);
        for (int j = 0; j < d; ++j) {
            outputFile << " " << A(i, j);
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

    bool Hpoly = false, Vpoly = false, Zono = false, cube = false, cross = false, simplex = false, prod_simplex = false,
            skinny_cube = false, randv = false;
    int d = 0, m = 0, kind = -1;

    for (int i = 1; i < argc; ++i) {
        bool correct = false;

        if (!strcmp(argv[i], "-help")) {
            std::cerr <<
                      "Usage:\n" <<
                      "-rz : generate a random zonotope\n" <<
                      "-dist : distribution to pick the length from [0,100] of each segment of a random zonotope: '1' for uniform , '2' for gaussian, '3' for exponential. The default generator is the uniform\n"
                      <<
                      "-rvc : generate a random V-polytope. The generator Samples uniformly distributed vertices from the boundary of a hypershere\n"
                      <<
                      "-rvs : generate a random V-polytope. The generator Samples uniformly distributed vertices from the interior of the unit cube\n"
                      <<
                      "-rhs : generate a random H-polytope\n" <<
                      "-cube : generate a hypercube\n" <<
                      "-cross : generate a cross polytope\n" <<
                      "-simplex : generate a simplex\n" <<
                      "-prod_simplex : generate a product of two simplices\n" <<
                      "-skinny_cube : generate a skinny hypercube\n" <<
                      "-h : generate a known polytope in H-representation\n" <<
                      "-v : generate a known polytope in V-representation\n" <<
                      "-d : the dimension\n" <<
                      "-m : number of zonotope's segments or V-polytope's vertices or H-polytope's facets\n" <<
                      std::endl;
            return 0;
        }

        if (!strcmp(argv[i], "-zonotope")) {
            Zono = true;
            kind = 1;
            correct = true;
        }
        if (!strcmp(argv[i], "-cube")) {
            cube = true;
            kind = 1;
            correct = true;
        }
        if (!strcmp(argv[i], "-cross")) {
            cross = true;
            kind = 2;
            correct = true;
        }
        if (!strcmp(argv[i], "-simplex")) {
            simplex = true;
            kind = 3;
            correct = true;
        }
        if (!strcmp(argv[i], "-prod_simplex")) {
            prod_simplex = true;
            kind = 4;
            correct = true;
        }
        if (!strcmp(argv[i], "-skinny_cube")) {
            skinny_cube = true;
            kind = 5;
            correct = true;
        }
        if (!strcmp(argv[i], "-h")) {
            Hpoly = true;
            correct = true;
        }
        if (!strcmp(argv[i], "-v")) {
            Vpoly = true;
            correct = true;
        }
        if (!strcmp(argv[i], "-rvc")) {
            Vpoly = true;
            kind = 4;
            randv = true;
            correct = true;
        }
        if (!strcmp(argv[i], "-rvs")) {
            Vpoly = true;
            kind = 5;
            randv = true;
            correct = true;
        }
        if (!strcmp(argv[i], "-rz")) {
            Zono = true;
            correct = true;
        }
        if (!strcmp(argv[i], "-rhs")) {
            Hpoly = true;
            kind = 6;
            correct = true;
        }
        if (!strcmp(argv[i], "-d")) {
            d = atof(argv[++i]);
            correct = true;
        }
        if (!strcmp(argv[i], "-dist")) {
            kind = atof(argv[++i]);
            correct = true;
        }
        if (!strcmp(argv[i], "-m")) {
            m = atof(argv[++i]);
            correct = true;
        }

        if (correct == false) {
            std::cerr << "unknown parameters \'" << argv[i] <<
                      "\', try " << argv[0] << " -help" << std::endl;
            exit(-2);
        }

    }

    if (d < 1) {
        std::cout << "Wrong inputs, try -help" << std::endl;
        exit(-1);
    }

    if (Zono) {
        if (kind > 3 || kind <1) std::cout << "Wrong declaration of zonotope's generator, try -help" << std::endl;
        if (m > 0) {
            Zonotope ZP;
            switch (kind) {
                case 1:
                    ZP = gen_zonotope_uniform<Zonotope, RNGType>(d, m);
                    break;
                case 2:
                    ZP = gen_zonotope_gaussian<Zonotope, RNGType>(d, m);
                    break;
                case 3:
                    ZP = gen_zonotope_exponential<Zonotope, RNGType>(d, m);
                    break;
                default:
                    kind = 1;
                    ZP = gen_zonotope_uniform<Zonotope, RNGType>(d, m);
                    break;
            }
            create_txt(ZP.get_mat(), ZP.get_vec(), kind, false, true);
        } else {
            std::cout << "Wrong inputs, try -help" << std::endl;
            exit(-1);
        }
    } else if (Hpoly) {
        Hpolytope HP;
        switch (kind) {
            case 1:
                HP = gen_cube<Hpolytope>(d, false);
                break;
            case 2:
                HP = gen_cross<Hpolytope>(d, false);
                break;
            case 3:
                HP = gen_simplex<Hpolytope>(d, false);
                break;
            case 4:
                HP = gen_prod_simplex<Hpolytope>(d);
                break;
            case 5:
                HP = gen_skinny_cube<Hpolytope>(d);
                break;
            case 6:
                if (m > d + 1) {
                    HP = random_hpoly<Hpolytope, RNGType>(d, m);
                } else {
                    std::cout << "The number of facets has to be >= d+1" << std::endl;
                    exit(-1);
                }
                break;
            default:
                std::cout << "Wrong inputs, try -help" << std::endl;
                exit(-1);
        }
        create_txt(HP.get_mat(), HP.get_vec(), kind, false, false);
    } else if (Vpoly) {
        Vpolytope VP;
        switch (kind) {
            case 1:
                VP = gen_cube<Vpolytope>(d, true);
                break;
            case 2:
                VP = gen_cross<Vpolytope>(d, true);
                break;
            case 3:
                VP = gen_simplex<Vpolytope>(d, true);
                break;
            case 4:
                if (!randv) {
                    std::cout<<"No prod_simplex in V-representation can be generated, try -help"<<std::endl;
                    exit(-1);
                }
                if (m > d) {
                    VP = random_vpoly_incube<Vpolytope, RNGType>(d, m);
                } else {
                    std::cout << "The number of vertices has to be >=d+1" << std::endl;
                    exit(-1);
                }
                break;
            case 5:
                if (!randv) {
                    std::cout<<"No skinny_cube in V-representation can be generated, try -help"<<std::endl;
                    exit(-1);
                }
                if (m > d) {
                    VP = random_vpoly<Vpolytope, RNGType>(d, m);
                } else {
                    std::cout << "The number of vertices has to be >=d+1" << std::endl;
                    exit(-1);
                }
                break;
            default:
                std::cout << "Wrong inputs, try -help" << std::endl;
                exit(-1);
        }
        create_txt(VP.get_mat(), VP.get_vec(), kind, true, false);
    } else {
        std::cout << "Wrong inputs, try -help" << std::endl;
        exit(-1);
    }
    return 0;
}
