#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "hpolytope.h"
#include "known_polytope_generators.h"

#include "random_walks/random_walks.hpp"

#include "volume_sequence_of_balls.hpp"
#include "volume_cooling_gaussians.hpp"
#include "volume_cooling_balls.hpp"

#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include "misc.h"


enum VOL_OPTIONS {
    SOB,
    CG,
    CB
};

enum ROUND_OPTIONS {
    SVD,
    MIN_ELLIPSOID
};


typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
typedef HPolytope <Point> HPOLYTOPE;
typedef typename HPOLYTOPE::MT MT;
typedef typename HPOLYTOPE::VT VT;

struct ArgOptions {
    VOL_OPTIONS vo;
    ROUND_OPTIONS ro;
    bool with_rounding;
    HPOLYTOPE* HP;

    NT volume_method (HPOLYTOPE& P, 
                      NT e, 
                      const unsigned int& walk_len) {
        switch (vo) {
            case SOB:
                return volume_sequence_of_balls<CDHRWalk, RNGType>(P, e, walk_len);;
            case CG:
                return volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(P, e, walk_len);
            case CB:
                return volume_cooling_balls<CDHRWalk, RNGType>(P, e, 2*walk_len).second;;
        }

        return -1;
    }

    std::tuple<MT, VT, NT> rounding_method (HPOLYTOPE& P, 
                                            std::pair<Point, NT>& InnerBall, 
                                            const unsigned int& walk_len, 
                                            RNGType& rng) {
        switch (ro) {
            case SVD:
                return svd_rounding<CDHRWalk, MT, VT>(P, InnerBall, walk_len, rng);
            case MIN_ELLIPSOID:
                return min_sampling_covering_ellipsoid_rounding<CDHRWalk, MT, VT>(P, InnerBall, walk_len, rng);
        }
    }
};

NT calculateLinearExtension(ArgOptions& args) {
    // Setup parameters for calculating volume and rounding
    unsigned int d = (args.HP)->dimension();
    unsigned int walk_len = 10 + d/10;
    NT e=0.1;
    // calculate volume of the order polytope
    NT round_multiply = 1.0;
    if (args.with_rounding) {
        //walk_len = 1;

        RNGType rng(d);
        std::pair<Point, NT> InnerBall = (args.HP)->ComputeInnerBall();
        std::tuple<MT, VT, NT> res = args.rounding_method(*(args.HP), InnerBall, 10 + 10*d, rng);
        round_multiply = std::get<2>(res);
    }
    NT volume = round_multiply * args.volume_method(*(args.HP), e, walk_len);

    // multiplying by d factorial, d = number of elements
    for(NT i=(NT)d; i>1; i-=1) {
        volume = volume * i;
    }
    return volume;
}


bool parseArgs(int argc, char* argv[], ArgOptions& args) {
    if(argc < 3) {
        std::cerr << "Too few arguments";
        std::cerr << "Usage: ./volesti_lecount INSTANCE VOLUME_METHOD ROUNDING_METHOD(optional) ";
        return false;
    }

    // create volume method
    std::string vm (argv[2]);
    if (vm == "cb") {
        args.vo = CB;
    }
    else if (vm == "cg") {
        args.vo = CG;
    }
    else if (vm == "sob") {
        args.vo = SOB;
    }
    else {
        std::cerr << "Invalid option for volume method";
        return false;
    }

    // create rounding method
    if (argc == 4) {
        std::string rm (argv[3]);
        if (rm == "SVD") {
            args.ro = SVD;
        }
        else if(rm == "MIN_ELLIPSOID") {
            args.ro = MIN_ELLIPSOID;
        }
        else {
            std::cerr << "Invalid option for rounding method";
            return false;
        }
        args.with_rounding = true;
    }
    else {
        args.with_rounding = false;    
    }

    // ----- START: parse the instance file and create an order polytope ------
    std::string filename (argv[1]);
    std::ifstream in(filename);

    // for storing relations a<=b
    std::vector<std::pair<int, int>> edges;
    int n = 0;
    int x;

    // read a single line
    std::string line;
    std::getline(in, line);
    std::stringstream line_ss(line);
    while(line_ss >> x) {
        if(x) {
            edges.emplace_back(0, n);
        }
        ++n;
    }

    // read rest of the lines
    for(int a = 1; a < n; ++a) {
        for(int b = 0; b < n; ++b) {
            if(!(in >> x)) {
                std::cerr << "Invalid adjacency matrix";
                return false;
            }

            if(x) {
                edges.emplace_back(a, b);
            }
        }
    }

    MT A = Eigen::MatrixXd::Zero(2*n + edges.size(), n);
    VT b = Eigen::MatrixXd::Zero(2*n + edges.size(), 1);

    // first add (ai >= 0) or (-ai <= 0) rows    
    A.topLeftCorner(n, n) = -Eigen::MatrixXd::Identity(n, n);

    // next add (ai <= 1) rows
    A.block(n, 0, n, n) = Eigen::MatrixXd::Identity(n, n);
    b.block(n, 0, n, 1) = Eigen::MatrixXd::Constant(n, 1, 1.0);


    // next add the relations
    for(int idx=0; idx<edges.size(); ++idx) {
        std::pair<int, int> edge  = edges[idx];
        A(2*n + idx, edge.first)  = 1;
        A(2*n + idx, edge.second) = -1;
    }

    args.HP = new HPOLYTOPE(n, A, b);
//    (args.HP)->print();
    return true;
}


/**
 
 Usage: ./volesti_lecount INSTANCE VOLUME_METHOD ROUNDING_METHOD 
        (ROUNDING_METHOD is optional) 

 example: for (volume method = sequence of balls, rounding method = SVD)
    ./volesti_lecount instances/bipartite_0.5_008_0.txt sob SVD

*/
int main(int argc, char* argv[]) {
    ArgOptions args;
    if (!parseArgs(argc, argv, args)) {
        std::cout << "Errors occured";
        return 0;
    }       

    std::cout << calculateLinearExtension(args) << "\n";
    return 0;
}
