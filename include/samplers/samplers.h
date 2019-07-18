// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H

#include "Eigen"
#include "spectrahedron.h"

typedef double NT_MATRIX;
typedef Eigen::Matrix<NT_MATRIX,Eigen::Dynamic,Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT_MATRIX,Eigen::Dynamic,1> VT;


// Pick a random direction as a normilized vector
template <class RNGType, class Point, typename NT>
Point get_direction(unsigned int dim) {

    boost::normal_distribution<> rdist(0,1);
    std::vector<NT> Xs(dim,0);
    NT normal = NT(0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = rdist(rng);
        normal += Xs[i] * Xs[i];
    }
    normal=1.0/std::sqrt(normal);

    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = Xs[i] * normal;
    }
    Point p(dim, Xs.begin(), Xs.end());
    return p;
}


// Pick a random point from a d-sphere
template <class RNGType, class Point, typename NT>
Point get_point_on_Dsphere(unsigned int dim, NT radius){
    Point p = get_direction<RNGType, Point, NT>(dim);
    p = (radius == 0) ? p : radius * p;
    return p;
}


// Pick a random point from a d-ball
template <class RNGType, class Point, typename NT>
Point get_point_in_Dsphere(unsigned int dim, NT radius){

    boost::random::uniform_real_distribution<> urdist(0,1);
    NT U;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng2(seed);
    Point p = get_direction<RNGType, Point, NT>(dim);
    U = urdist(rng2);
    U = std::pow(U, 1.0/(NT(dim)));
    p = (radius*U)*p;
    return p;
}

// ball walk with uniform target distribution
template <class RNGType, class Point, class Polytope, typename NT>
void ball_walk(Point &p,
               Polytope &P,
               NT delta)
{
    //typedef typename Parameters::RNGType RNGType;
    Point y = get_point_in_Dsphere<RNGType, Point>(p.dimension(), delta);
    y = y + p;
    if (P.is_in(y)==-1) p = y;
}

// WARNING: USE ONLY WITH BIRKHOFF POLYOPES
// Compute more random points using symmetries of birkhoff polytope
/*
template <class T, class Point, class K, typename  NT>
int birk_sym(T &P, K &randPoints, Point &p) {
    int n=std::sqrt(p.dimension());
    std::vector<int> myints;
    for (unsigned int i=0; i<n; i++){
        myints.push_back(i);
    }

    //The n! possible permutations with n elements
    do {
        std::vector<NT> newpv;

        for (unsigned int j=0; j<p.dimension(); j++){
            int idx = (myints[j/n])*n+1+j%n-1;
            newpv.push_back(p[idx]);
        }

        Point new_p(p.dimension(),newpv.begin(),newpv.end());
        if(P.is_in(new_p) != 0) {
            randPoints.push_back(new_p);
        }
    } while ( std::next_permutation(myints.begin(),myints.end()) );
}*/


// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <class Polytope, class PointList, class Parameters, class Point>
void rand_point_generator(Polytope &P,
                         Point &p,   // a point to start
                         unsigned int rnum,
                         unsigned int walk_len,
                         PointList &randPoints,
                         Parameters &var)  // constants for volume
{
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    std::vector <NT> lamdas(P.num_of_hyperplanes(), NT(0));
    unsigned int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta;
    Point p_prev = p;

    if (var.ball_walk) {
        ball_walk <RNGType> (p, P, ball_rad);
    }else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else
        hit_and_run(p, P, var);

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk<RNGType>(p, P, ball_rad);
            }else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas);
            } else
                hit_and_run(p, P, var);
        }
        randPoints.push_back(p);
    }
 
}






template <class BallPoly, class PointList, class Parameters, class Point>
void rand_point_generator(BallPoly &PBLarge,
                         Point &p,   // a point to start
                         unsigned int rnum,
                         unsigned int walk_len,
                         PointList &randPoints,
                         BallPoly &PBSmall,
                         unsigned int &nump_PBSmall,
                         Parameters &var) {  // constants for volume

    typedef typename Point::FT NT;
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    std::vector <NT> lamdas(PBLarge.num_of_hyperplanes(), NT(0));
    unsigned int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta;
    Point p_prev = p;

    if (var.ball_walk) {
        ball_walk<RNGType>(p, PBLarge, ball_rad);
    }else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        
        std::pair <NT, NT> bpair = PBLarge.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else {
        hit_and_run(p, PBLarge, var);
    }

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk<RNGType>(p, PBLarge, ball_rad);
            }else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, PBLarge, rand_coord, rand_coord_prev, kapa, lamdas);
            } else
                hit_and_run(p, PBLarge, var);
        }
        if (PBSmall.second().is_in(p) == -1) {//is in
            randPoints.push_back(p);
            ++nump_PBSmall;
        }
    }
}

// ----- HIT AND RUN FUNCTIONS ------------ //

//hit-and-run with random directions and update
template <class Polytope, class Point, class Parameters>
void hit_and_run(Point &p,
                Polytope &P,
                Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    std::pair <NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
}

template <class Point, class Parameters>
void hit_and_run(Point& point,
        Spectrahedron &spectrahedron,
        Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, double>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(pointVT, lVT);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    Point b1 = (min_plus * l) + point;
    Point b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

template <class Point, class Parameters>
void hit_and_run(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 VT& a,
                 double b) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, double>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(pointVT, lVT, a, b);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    Point b1 = (min_plus * l) + point;
    Point b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

template <class Point, class Parameters>
void hit_and_run_sampled_covariance_matrix(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 VT& a,
                 double b,
                 MT& covarianceMatrix) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, double>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    lVT = covarianceMatrix * lVT;
    l = Point(lVT);
    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(pointVT, lVT, a, b);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    Point b1 = (min_plus * l) + point;
    Point b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

//returns at b1, b2 the intersection points of the ray with the polytope
template<class Polytope, class Point, class Parameters>
void hit_and_run(Point &p,
                 Polytope &P,
                 Parameters &var,
                 Point &b1,
                 Point &b2) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    typedef typename Point::FT NT;

    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    std::pair<NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    b1 = (min_plus * l) + p;
    b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
}

template <class Point, class Parameters>
void hit_and_run(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 Point &b1,
                 Point &b2) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, double>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(pointVT, lVT);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    b1 = (min_plus * l) + point;
    b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}


template <class Point, class Parameters>
void hit_and_run(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 Point &b1,
                 Point &b2,
                 VT& a,
                 double b) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, double>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(pointVT, lVT, a, b);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    b1 = (min_plus * l) + point;
    b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

template <class Point, class Parameters>
void hit_and_run_sampled_covariance_matrix(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 Point &b1,
                 Point &b2,
                 VT& a,
                 double b,
                 MT& covarianceMatrix) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, double>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    lVT = covarianceMatrix * lVT;
    l = Point(lVT);
    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(pointVT, lVT, a, b);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    b1 = (min_plus * l) + point;
    b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

//hit-and-run with random directions and update
// copied this function - need it to return the end points of the segment
template<class Polytope, class Point, class Parameters>
void implicit_isotropization_hit_and_run(Point &p,
                                         Polytope &P,
                                         Parameters &var,
                                         Point &b1,
                                         Point &b2,
                                         MT &isotropic) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;

    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_point_on_Dsphere<RNGType, Point, NT>(n, 1.0);//get_direction<RNGType, Point, NT>(n);
    MT v(n, 1);
    v = isotropic * l.getCoefficients();
//
    for (int i=0 ; i<l.dimension() ; i++) {
        l.set_coord(i, v(i, 0));
        i++;//TODO clean up
    }
//    l.add(isotropic*l.getCoefficients());

    std::pair<NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    b1 = (min_plus * l) + p;
    b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
}

//return at bpair the lambdas of the end points
//hit-and-run with orthogonal directions and update
template <class Polytope, class Point, typename NT>
void hit_and_run_coord_update(Point &p,
                              Point &p_prev,
                              Polytope &P,
                              unsigned int rand_coord,
                              unsigned int rand_coord_prev,
                              NT kapa,
                              std::vector<NT> &lamdas,
                              std::pair <NT, NT>& bpair) {
    bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
}

//return at bpair the lambdas of the end points and the endpoints
//hit-and-run with orthogonal directions and update
template <class Polytope, class Point, typename NT>
void hit_and_run_coord_update_isotropic(Point &p,
                              Point &p_prev,
                              Polytope &P,
                              unsigned int rand_coord,
                              unsigned int rand_coord_prev,
                              NT kapa,
                              std::vector<NT> &lamdas,
                              std::pair <NT, NT>& bpair,
                              Point& p1,
                              Point& p2,
                              VT& isotropic) {
    Point v = Point(isotropic);
    bpair = P.line_intersect(p,v);
//    bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    NT min_plus = bpair.first;
    NT max_minus = bpair.second;
    p1 = (min_plus * v) + p;
    p2 = (max_minus * v) + p;
    p = (kapa * p1);
    p = ((1 - kapa) * p2) + p;
//    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
}


//hit-and-run with orthogonal directions and update
template <class Polytope, class Point, typename NT>
void hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             Polytope &P,
                             unsigned int rand_coord,
                             unsigned int rand_coord_prev,
                             NT kapa,
                             std::vector<NT> &lamdas) {
    std::pair <NT, NT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
}


#endif //RANDOM_SAMPLERS_H
