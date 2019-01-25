// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "polytopes.h"
#include "samplers.h"
#include "gaussian_samplers.h"
#include "sample_only.h"
#include "simplex_samplers.h"
#include "vpolyintersectvpoly.h"

/*
Rcpp::NumericMatrix Rsample_points (Rcpp::NumericMatrix A, unsigned int walk_len, Rcpp::NumericVector InnerPoint,
                                    bool gaussian, bool ball_walk, double delta, bool coord, bool Vpoly, bool Zono, bool sam_simplex,
                                    bool sam_can_simplex, bool sam_arb_simplex, bool sam_ball, bool sam_sphere,
                                    unsigned int numpoints, int dim, double variance) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> zonotope;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;
    InterVP VPcVP;

    std::list<Point> randPoints;
    Rcpp::NumericMatrix PointSet(dim,numpoints);

    if (sam_ball || sam_sphere) {

        for (unsigned int k = 0; k < numpoints; ++k) {
            if (sam_ball) {
                randPoints.push_back(get_point_in_Dsphere<RNGType , Point >(dim, delta));
            } else {
                randPoints.push_back(get_point_on_Dsphere<RNGType , Point >(dim, delta));
            }
        }

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        unsigned int j = 0, i;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;

    }

    if (sam_simplex || sam_can_simplex) {

        if (sam_simplex) {
            Sam_Unit<NT, RNGType >(dim, numpoints, randPoints);
        } else {
            Sam_Canon_Unit<NT, RNGType >(dim, numpoints, randPoints);
        }

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        unsigned int j = 0, i;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    if (sam_arb_simplex) {
        unsigned int n=A.ncol()-1;
        std::vector<NT> temp_p(n, 0.0);
        typename std::vector<NT>::iterator temp_it;
        std::vector<Point> vec_point;

        for (int k = 1; k < A.nrow(); ++k) {
            temp_it = temp_p.begin();
            for (int l = 1; l < A.ncol(); ++l, ++temp_it) {
                *temp_it = A(k,l);
            }
            vec_point.push_back(Point(n, temp_p.begin(), temp_p.end()));
        }

        Sam_arb_simplex<NT, RNGType>(vec_point.begin(), vec_point.end(), numpoints, randPoints);

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        unsigned int j = 0, i;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose=false,
            coordinate=coord;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    unsigned int m=A.nrow()-1;
    unsigned int n=A.ncol()-1;
    std::vector <std::vector<NT> > Pin(m + 1, std::vector<NT>(n + 1));
    Rcpp::NumericMatrix PointSet2(n,numpoints);

    for (unsigned int i = 0; i < m + 1; i++) {
        for (unsigned int j = 0; j < n + 1; j++) {
            Pin[i][j] = A(i, j);
        }
    }
    // construct polytope
    if (Zono) {
        ZP.init(Pin);
    } else if (!Vpoly) {
        HP.init(Pin);
    } else {
        VP.init(Pin);
    }

    std::pair<Point,NT> InnerBall;
    if (Zono) {
        InnerBall = ZP.ComputeInnerBall();
    } else if (!Vpoly) {
        InnerBall = HP.ComputeInnerBall();
    } else {
        InnerBall = VP.ComputeInnerBall();
    }

    if (InnerPoint.size()==n) {
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerPoint[j]);
        }
        InnerBall.first = Point( n , temp_p.begin() , temp_p.end() );
    }

    Point p = InnerBall.first;
    NT a = 1.0 / (2.0 * variance);
    if (ball_walk) {
        if (delta < 0.0) { // set the radius for the ball walk if is not set by the user
            if (gaussian) {
                delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
            } else {
                delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
            }
        }
    }
    unsigned int rnum = std::pow(1.0,-2) * 400 * n * std::log(n);
    // initialization
    vars<NT, RNGType> var1(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                           delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
    vars_g<NT, RNGType> var2(n, walk_len, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                             rand_only, false, NN, birk, ball_walk, coord);
    if (Zono) {
        sampling_only<Point>(randPoints, ZP, walk_len, numpoints, gaussian, a, p, var1, var2);
    } else if (!Vpoly) {
        sampling_only<Point>(randPoints, HP, walk_len, numpoints, gaussian, a, p, var1, var2);
    } else {
        sampling_only<Point>(randPoints, VP, walk_len, numpoints, gaussian, a, p, var1, var2);
    }

    // store the sampled points to the output matrixRcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["radius"])
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    unsigned int j = 0, i;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet2(i,j)=*qit;
        }
    }
    return PointSet2;

}*/


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_points(Rcpp::Nullable<Rcpp::Reference> P = R_NilValue, Rcpp::Nullable<unsigned int> N = R_NilValue,
                                  Rcpp::Nullable<std::string> WalkType = R_NilValue,
                                  Rcpp::Nullable<unsigned int> walk_len = R_NilValue,
                                  Rcpp::Nullable<bool> exact = R_NilValue,
                                  Rcpp::Nullable<std::string> body = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue,
                                  Rcpp::Nullable<std::string> distribution = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> InnerPoint = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;
    InterVP VPcVP;

    int type, dim, numpoints;
    NT radius = 1.0, delta = -1.0;
    bool set_mean_point = false, coordinate, ball_walk, gaussian = false;
    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall;

    if (!N.isNotNull()) {
        numpoints = 100;
    } else {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    if (exact.isNotNull()) {
        if (P.isNotNull()) {
            type = Rcpp::as<Rcpp::Reference>(P).field("t");
            if (Rcpp::as<bool>(exact) && type==2) {
                if (Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows() == Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols()+1) {
                    Vpolytope VP;
                    VP.init(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols(),
                            Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")),
                            VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows()));
                    Sam_arb_simplex(VP, numpoints, randPoints);
                } else {
                    throw Rcpp::exception("Not a simplex!");
                }
            } else if (Rcpp::as<bool>(exact) && type!=2) {
                throw Rcpp::exception("Not a simplex in V-representation!");
            }
        } else {
            if (!body.isNotNull()) {

                throw Rcpp::exception("Wrong input!");

            } else if (!Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("dimension")) {

                throw Rcpp::exception("Wrong input!");

            }
            dim = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["dimension"]);
            if (!Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("radius")) {

                radius = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["radius"]);

            }
            if (Rcpp::as<std::string>(body).compare(std::string("hypersphere"))==0) {

                for (unsigned int k = 0; k < numpoints; ++k) {
                    randPoints.push_back(get_point_on_Dsphere<RNGType , Point >(dim, radius));
                }

            } else if (Rcpp::as<std::string>(body).compare(std::string("ball"))==0) {

                for (unsigned int k = 0; k < numpoints; ++k) {
                    randPoints.push_back(get_point_in_Dsphere<RNGType , Point >(dim, radius));
                }

            } else if (Rcpp::as<std::string>(body).compare(std::string("unit simplex"))==0) {

                Sam_Unit<NT, RNGType >(dim, numpoints, randPoints);

            } else if (Rcpp::as<std::string>(body).compare(std::string("canonical simplex"))==0) {

                Sam_Canon_Unit<NT, RNGType >(dim, numpoints, randPoints);

            } else {

                throw Rcpp::exception("Wrong input!");

            }
        }
    } else if (P.isNotNull()) {

        type = Rcpp::as<Rcpp::Reference>(P).field("type");
        dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");
        unsigned int walkL = 10+dim/10;
        Point MeanPoint;
        if (InnerPoint.isNotNull()) {
            if (Rcpp::as<Rcpp::NumericVector>(InnerPoint).size()!=dim) {
                Rcpp::warning("Internal Point has to lie in the same dimension as the polytope P");
            } else {
                set_mean_point = true;
                MeanPoint = Point( dim , Rcpp::as<std::vector<NT> >(InnerPoint).begin(),
                        Rcpp::as<std::vector<NT> >(InnerPoint).end() );
            }
        }
        if(walk_len.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_len);

        NT a = 0.5;

        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("variance"))
            a = 1.0 / (2.0 * Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["variance"]));

        if(!WalkType.isNotNull() || Rcpp::as<std::string>(WalkType).compare(std::string("CDHR"))==0){
            coordinate = true;
            ball_walk = false;
        } else if (Rcpp::as<std::string>(WalkType).compare(std::string("RDHR"))==0) {
            coordinate = false;
            ball_walk = false;
        } else if (Rcpp::as<std::string>(WalkType).compare(std::string("BW"))==0) {
            if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("BW_rad")) {
                delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["BW_rad"]);
            }
            coordinate = false;
            ball_walk = true;
        } else {
            throw Rcpp::exception("Unknown walk type!");
        }

        if (distribution.isNotNull()) {
            if (Rcpp::as<std::string>(distribution).compare(std::string("gaussian"))==0) {
                gaussian = true;
            } else if(Rcpp::as<std::string>(distribution).compare(std::string("uniform"))!=0) {
                throw Rcpp::exception("Wrong distribution!");
            }
        }
        bool rand_only=false,
                NN=false,
                birk=false,
                verbose=false;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // the random engine with this seed
        RNGType rng(seed);
        boost::random::uniform_real_distribution<>(urdist);
        boost::random::uniform_real_distribution<> urdist1(-1,1);
        vars<NT, RNGType> var1(1,dim,walkL,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                               delta,verbose,rand_only,false,NN,birk,ball_walk,coordinate);
        vars_g<NT, RNGType> var2(dim, walkL, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                                 rand_only, false, NN, birk, ball_walk, coordinate);

        if (type==1) {
            // Hpolytope
            //Hpolytope HP;
            HP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("A")), Rcpp::as<VT>(Rcpp::as<Rcpp::Reference>(P).field("b")));

            if (!set_mean_point || ball_walk) {
                InnerBall = HP.ComputeInnerBall();
                if (!set_mean_point) MeanPoint = InnerBall.first;
            }

        } else if(type==2) {
            // Vpolytope
            //Vpolytope VP;
            VP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")), VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows()));

            if (!set_mean_point || ball_walk) {
                InnerBall = VP.ComputeInnerBall();
                if (!set_mean_point) MeanPoint = InnerBall.first;
            }

        } else if(type==3){
            // Zonotope
            //zonotope ZP;
            ZP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")), VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));

            if (!set_mean_point || ball_walk) {
                InnerBall = VP.ComputeInnerBall();
                if (!set_mean_point) MeanPoint = InnerBall.first;
            }

        } else {
            // Intersection of two V-polytopes
            Vpolytope VP1;
            Vpolytope VP2;
            //InterVP VPcVP;
            VP1.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V1")), VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V1")).rows()));
            VP2.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V2")), VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V2")).rows()));
            VPcVP.init(VP1, VP2);

            if (!set_mean_point || ball_walk) {
                InnerBall = VP.ComputeInnerBall();
                if (!set_mean_point) MeanPoint = InnerBall.first;
            }
        }

        if (ball_walk) {
            if (gaussian) {
                delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(dim));
            } else {
                delta = 4.0 * InnerBall.second / std::sqrt(NT(dim));
            }
        }

        if (type == 1) {
            sampling_only<Point>(randPoints, HP, walkL, numpoints, gaussian,
                                 a, MeanPoint, var1, var2);
        } else if (type == 2) {
            sampling_only<Point>(randPoints, VP, walkL, numpoints, gaussian,
                                 a, MeanPoint, var1, var2);
        } else if (type == 3) {
            sampling_only<Point>(randPoints, ZP, walkL, numpoints, gaussian,
                                 a, MeanPoint, var1, var2);
        } else {
            sampling_only<Point>(randPoints, VPcVP, walkL, numpoints, gaussian,
                                 a, MeanPoint, var1, var2);
        }

    } else {

        throw Rcpp::exception("Wrong input!");

    }

    Rcpp::NumericMatrix PointSet(dim,numpoints);
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    unsigned int j = 0, i;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet(i,j)=*qit;
        }
    }
    return PointSet;

}
