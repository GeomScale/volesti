// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
#include <Rcpp.h>
#include <RcppEigen.h>
#include "volume.h"
#include "extractMatPoly.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix vol_R (Rcpp::NumericMatrix A, int walk_len, double e, Rcpp::NumericVector Chebychev,
                           bool annealing, int win_len, int N, double C, double ratio, double frac,
                           bool ball_walk, double delta, bool Vpoly, bool Zono, bool exact_zono, bool gen_only,
                           bool Vpoly_gen, int kind_gen, int dim_gen, int m_gen, bool round_only,
                           bool rotate_only, bool ball_only, bool sample_only, int numpoints, double variance,
                           bool coord, bool rounding, bool verbose) {


    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> Zonotope;
    int nexp=1, n_threads=1,i,j;
    //NT exactvol(-1.0);
    bool rand_only=false,
	 file=false,
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         experiments=true,
         coordinate=coord;
    Hpolytope HP;
    Vpolytope VP;
    Zonotope ZP;

    int m=A.nrow()-1;
    int n=A.ncol()-1;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    Rcpp::NumericMatrix vol_res(1,1);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    //boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));
    std::vector<NT> bin(m);

    if (gen_only) {
        Rcpp::NumericMatrix Mat;
        if (kind_gen == 0) {
            Zonotope ZP = gen_zonotope<Zonotope, RNGType>(dim_gen, m_gen);
            Mat = extractMatPoly(ZP);
        } else if (Vpoly_gen) {
            if (kind_gen == 1) {
                VP = gen_cube<Vpolytope>(dim_gen, true);
                Mat = extractMatPoly(VP);
            } else if (kind_gen == 2) {
                VP = gen_cross<Vpolytope>(dim_gen, true);
                Mat = extractMatPoly(VP);
            } else if (kind_gen == 3) {
                VP = gen_simplex<Vpolytope>(dim_gen, true);
                Mat = extractMatPoly(VP);
            } else {
                std::cout<<"An internal error is occured.. We are sorry for the inconvinience."<<std::endl;
                exit(-1);
            }
        } else {
            if (kind_gen == 1) {
                HP = gen_cube<Hpolytope>(dim_gen, false);
                Mat = extractMatPoly(HP);
            } else if (kind_gen == 2) {
                HP = gen_cross<Hpolytope>(dim_gen, false);
                Mat = extractMatPoly(HP);
            } else if (kind_gen == 3) {
                HP = gen_simplex<Hpolytope>(dim_gen, false);
                Mat = extractMatPoly(HP);
            } else if (kind_gen == 4) {
                HP = gen_prod_simplex<Hpolytope>(dim_gen);
                Mat = extractMatPoly(HP);
            } else if (kind_gen == 5) {
                HP = gen_skinny_cube<Hpolytope>(dim_gen);
                Mat = extractMatPoly(HP);
            } else {
                std::cout<<"An internal error is occured.. We are sorry for the inconvinience."<<std::endl;
                exit(-1);
            }
        }
        return Mat;
    }

    for (i=0; i<m+1; i++){
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }
    // construct polytope
    if (Zono) {
        ZP.init(Pin);
        if (exact_zono) {
            vol_res(0,0) = exact_zonotope_vol<NT>(ZP);
            return vol_res;
        }
    } else if (!Vpoly) {
        HP.init(Pin);
    } else {
        VP.init(Pin);
    }

    std::pair<Point,NT> InnerBall;
    if (ball_only) {
        Rcpp::NumericMatrix Mat(1, n + 1);
        if (Zono) {
            InnerBall = ZP.ComputeInnerBall();
        } else if (!Vpoly) {
            InnerBall = HP.ComputeInnerBall();
        } else {
            InnerBall = VP.ComputeInnerBall();
        }
        for (int k = 0; k < n; ++k) {
            Mat(0,k) = InnerBall.first[k];
        }
        Mat(0,n) = InnerBall.second;
        return Mat;
    }

    if (rotate_only) {
        Rcpp::NumericMatrix Mat;
        if (Zono) {
            rotating<NT>(ZP);
            Mat = extractMatPoly(ZP);
        }else if (!Vpoly) {
            rotating<NT>(HP);
            Mat = extractMatPoly(HP);
        } else {
            rotating<NT>(VP);
            Mat = extractMatPoly(VP);
        }
        return Mat;
    }


    //Compute chebychev ball//
    if(Chebychev.size()==n+1 || Chebychev.size()==n) { //if it is given as an input
        if (Chebychev.size()==n) {
            // if only sampling is requested
            // the radius of the inscribed ball is going to be needed for the sampling (radius of ball walk)
            if (Zono) {
                InnerBall = ZP.ComputeInnerBall();
            } else if (!Vpoly) {
                InnerBall = HP.ComputeInnerBall();
            } else {
                InnerBall = VP.ComputeInnerBall();
            }
        }
        // store internal point hat is given as input
        std::vector<NT> temp_p;
        for (int j=0; j<n; j++){
            temp_p.push_back(Chebychev[j]);
        }
        InnerBall.first = Point( n , temp_p.begin() , temp_p.end() );
        // store the radius of the internal ball that is given as input
        if (Chebychev.size()==n+1) InnerBall.second = Chebychev[n];
    } else {
        // no internal ball or point is given as input
        if (Zono) {
            InnerBall = ZP.ComputeInnerBall();
        } else if (!Vpoly) {
            InnerBall = HP.ComputeInnerBall();
        } else {
            InnerBall = VP.ComputeInnerBall();
        }
    }

    // if only rounding is requested
    if (round_only) {
        Rcpp::NumericMatrix Mat;
        if (ball_walk) {
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
        }
        // initialization
        vars<NT, RNGType> var(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                 delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
        std::pair <NT, NT> round_res;
        if (Zono) {
            round_res = rounding_min_ellipsoid(ZP, InnerBall, var);
            Mat = extractMatPoly(ZP);
        } else if (!Vpoly) {
            round_res = rounding_min_ellipsoid(HP, InnerBall, var);
            Mat = extractMatPoly(HP);
        } else {
            round_res = rounding_min_ellipsoid(VP, InnerBall, var);
            Mat = extractMatPoly(VP);
        }
        // store rounding value and the ratio between min and max axe in the first row
        // the matrix is in ine format so the first row is useless and is going to be removed by R function modifyMat()
        Mat(0,0) = round_res.first;
        Mat(0,1) = round_res.second;
        return Mat;
    }

    // if only sampling is requested
    if (sample_only){
        std::list<Point> randPoints;
        Point p = InnerBall.first;
        NT a = 1.0 / (2.0 * variance);
        if (ball_walk){
            if(delta<0.0){ // set the radius for the ball walk if is not set by the user
                if(annealing) {
                    delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
                } else {
                    delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
                }
            }
        }
        // initialization
        vars<NT, RNGType> var1(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                 delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
        vars_g<NT, RNGType> var2(n, walk_len, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                    rand_only, false, NN, birk, ball_walk, coord);
        if (Zono) {
            sampling_only<Point>(randPoints, ZP, walk_len, numpoints, annealing, a, p, var1, var2);
        } else if(!Vpoly) {
            sampling_only<Point>(randPoints, HP, walk_len, numpoints, annealing, a, p, var1, var2);
        } else {
            sampling_only<Point>(randPoints, VP, walk_len, numpoints, annealing, a, p, var1, var2);
        }
        Rcpp::NumericMatrix PointSet(n,numpoints);

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        j = 0;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    // print chebychev ball in verbose mode
    if (verbose) {
        std::cout << "Inner ball center = " << std::endl;
        for (i = 0; i < n; i++) {
            std::cout << InnerBall.first[i] << " ";
        }
        std::cout << "\nradius of inner ball = " << InnerBall.second << std::endl;
    }


    // initialization
    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0, InnerBall.second,rng,urdist,urdist1,
             delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (annealing) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                  urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g<NT, RNGType> var1(n, walk_len, N, win_len, 1, e, InnerBall.second, rng, C, frac, ratio, delta, false, verbose,
                    rand_only, rounding, NN, birk, ball_walk, coordinate);
        if (Zono) {
            vol = volume_gaussian_annealing(ZP, var1, var2, InnerBall);
        } else if (!Vpoly) { // if the input is a H-polytope
            vol = volume_gaussian_annealing(HP, var1, var2, InnerBall);
        } else {  // if the input is a V-polytope
            vol = volume_gaussian_annealing(VP, var1, var2, InnerBall);
        }
        if (verbose) std::cout << "volume computed = " << vol << std::endl;
    } else {
        if (Zono) {
            vol = volume(ZP, var, var, InnerBall);
        } else if (!Vpoly) { // if the input is a H-polytope
            vol = volume(HP, var, var, InnerBall);
        } else { // if the input is a V-polytope
            vol = volume(VP, var, var, InnerBall);
        }
    }
    vol_res(0,0) = vol;
    
    return vol_res;
}
