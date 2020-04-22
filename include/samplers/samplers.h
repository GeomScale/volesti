// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H

//#include "Eigen"
//#include "spectrahedron.h"
#include <cmath>

//extern int BOUNDARY_CALLS;
//extern double ORACLE_TIME;
//extern double REFLECTION_TIME;

typedef double NT_MATRIX;
//typedef Eigen::Matrix<NT_MATRIX,Eigen::Dynamic,Eigen::Dynamic> MT;
//typedef Eigen::Matrix<NT_MATRIX,Eigen::Dynamic,1> VT;


// Pick a random direction as a normilized vector
template <class RNGType, class Point, typename NT>
Point get_direction(unsigned int dim) {

    boost::normal_distribution<> rdist(0,1);
    std::vector<NT> Xs(dim,0);
    NT normal = NT(0);
    unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();//4 if fixed for debug
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




/**
 * Return a random vector v according to the n-dimensional normal distribution with mean 0 and covariance
 * matrix V
 *
 * @tparam RNGType
 * @tparam Point
 * @tparam NT
 * @param dim
 * @param choleskyDecomp the cholesky decomposition of V
 * @return
 */
template <class MT, class RNGType, class Point, typename NT>
Point get_direction(unsigned int dim, const MT& choleskyDecomp) {
    Point l = get_direction<RNGType, Point, NT>(dim);
    return Point(choleskyDecomp * l.getCoefficients());
}

// Pick a random point from a d-sphere
template <class RNGType, class Point, typename NT>
Point get_point_on_Dsphere(unsigned int dim, NT radius){
    Point p = get_direction<RNGType, Point, NT>(dim);
    p = (radius == 0) ? p : radius * p;
    return p;
}

/**
 * direction for "running Shake-and-Bake" as in 1.3.3 of
 * Boender et al. (1991)
 */
template <class RNGType, class Point, typename NT, class Parameters>
void rsab_get_direction(const unsigned int& dim, const Point& normal, Point& direction, Parameters& var) {
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);


    NT r = pow(urdist(rng), 1.0/(dim - 1));
//    NT r = pow(0.5, 1.0/(dim - 1));
    direction = get_point_on_Dsphere<RNGType, Point, NT> (dim,1);

    NT cd = direction.dot(normal);
    NT fd = r / sqrt(1 - cd * cd);
    NT fc = -(r * cd / sqrt(1 - cd * cd) + sqrt(1 - r * r));
    direction *= fd ;
    direction += fc * normal;
    direction.normalize();
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
                         Parameters const& var)  // constants for volume
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

template <class Spectrahedron, class Parameters, class SpecSettings, class Point>
void uniform_next_point_spec(Spectrahedron &spectrahedron,
                          Point &p,   // a point to start
                          unsigned int walk_len,
                          Parameters &var,
                          SpecSettings &settings) {

    var.steps = var.steps + walk_len;
    for (unsigned int j = 0; j < walk_len; ++j) {
        if (var.billiard) {
            billiard_walk(spectrahedron, p, var.diameter, var, p, 0.0,
                             settings, false);
        } else {
            //hit_and_run_spec(p, spectrahedron, var);
        }
    }
}


template <class Spectrahedron, class PointList, class Parameters, class SpecSettings, class Point>
void rand_point_generator_spec(Spectrahedron &spectrahedron,
                          Point &p,   // a point to start
                          unsigned int rnum,
                          unsigned int walk_len,
                          PointList &randPoints,
                          Parameters &var,
                          SpecSettings &settings)  {

    var.steps = var.steps + rnum * walk_len;

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.billiard) {
                billiard_walk(spectrahedron, p, var.diameter, var, p, 0.0,
                                 settings, false);
            } else {
                //hit_and_run_spec(p, spectrahedron, var);
            }
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
                         Parameters const& var) {  // constants for volume

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


template<class Polytope, class Parameters, class Point, class MT>
void rand_point_generator_Boltzmann(Polytope &P,
                                        Point& c,
                                        Point &p,   // a point to start
                                        unsigned int walk_length,
                                        Parameters &var,
                                        double temperature,
                                        const MT& covariance_matrix) {

    // begin sampling
    for (unsigned int i = 1; i <= walk_length ; ++i) {
        hit_and_run_Boltzmann(p, P, var, c, temperature, covariance_matrix);
    }
}

template<class Polytope, class Parameters, class Point, class MT, class PointsList>
void rand_point_generator_Boltzmann(Polytope &P,
                                    Point& c,
                                    Point &p,   // a point to start
                                    unsigned int pointsNum,
                                    unsigned int walk_length,
                                    Parameters &var,
                                    double temperature,
                                    const MT& covariance_matrix,
                                    PointsList& points) {

    // begin sampling
    Point temp = p;

    for (int at=0 ; at<pointsNum ; at++) {
        p = temp;

        for (unsigned int i = 1; i <= walk_length; ++i) {
            hit_and_run_Boltzmann(p, P, var, c, temperature, covariance_matrix);
        }

        points.push_back(p);
    }
}



template <class Spectrahedron, class PointList, class Parameters, class Point>
void boundary_rand_point_generator_spec(Spectrahedron &spectrahedron,
                                   Point &p,   // a point to start
                                   unsigned int rnum,
                                   unsigned int walk_len,
                                   PointList &randPoints,
                                   Parameters &var)  // constants for volume
{
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    typedef typename Spectrahedron::VT VT;
    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    //std::vector <NT> lamdas(P.num_of_hyperplanes(), NT(0)), Av(P.num_of_hyperplanes(), NT(0));
    // unsigned int rand_coord, rand_coord_prev;
    NT lambda;
    Point p1(n), p2(n), v(n);
    VT pointVT(n), lVT(n);
    std::pair <NT, NT> dbpair;


    //hit_and_run(p, P, var);

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {

            v = get_direction<RNGType, Point, NT>(n);
            pointVT = p.getCoefficients();
            lVT = v.getCoefficients();
            dbpair = spectrahedron.boundaryOracle(pointVT, lVT);
            p1 = (dbpair.first * v) + p;
            p2 = (dbpair.second * v) + p;
            lambda = urdist(rng) * (dbpair.first - dbpair.second) + dbpair.second;
            p = (lambda * v) + p;

        }
        randPoints.push_back(p1);
        randPoints.push_back(p2);
    }

}




// ----- HIT AND RUN FUNCTIONS ------------ //

//hit-and-run with random directions and update
template <class Polytope, class Point, class Parameters>
void hit_and_run(Point &p,
                Polytope &P,
                Parameters const& var) {
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

template <typename NT>
NT expCDF(NT x, NT lambda) {
  return 1.0 - std::exp(- lambda*x);
}

template <typename NT>
NT expQuantile(NT p, NT lambda) {
  return -std::log(1.0 - p) / lambda;
}

template <class RNGType, typename NT>
NT texp(NT lambda, NT a, NT b, RNGType& rng) {
  boost::random::uniform_real_distribution<> urdist(0, 1);
  NT u = urdist(rng);
  NT cdfA = expCDF(a, lambda);
  NT cdfB = expCDF(b, lambda);
  
  return expQuantile(cdfA + u*(cdfB - cdfA), lambda);
}

template <class Polytope, class Point, class Parameters, class MT>
void hit_and_run_Boltzmann(Point &p,
                           Polytope &P,
                           Parameters &var,
                           Point& BoltzmannDirection,
                           double BoltzmannParameter,
                           const MT& choleskyDecomp) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n =p.dimension();
    RNGType &rng = var.rng;

    Point l = get_direction<RNGType, Point, NT>(n, choleskyDecomp);

    std::pair <NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
    double c1 = BoltzmannDirection.dot(b1);
    double c2 = BoltzmannDirection.dot(b2);

    double lambda;

    if (c1 > c2) {
        lambda = texp((c1 - c2) / BoltzmannParameter, NT(0), min_plus - max_minus, rng);
        p = b2;
    }
    else {
        lambda = -texp((c2 - c1) / BoltzmannParameter, NT(0), min_plus - max_minus, rng);
        p = b1;
    }

    p = (lambda * l) + p;
}


template <class Spectrahedron, class Point, class Parameters, typename NT>
void hit_and_run_Boltzmann_spec(Point &p,
                                Spectrahedron &P,
                           Parameters &var,
                           Point& BoltzmannDirection,
                           NT BoltzmannParameter) {
    typedef typename Parameters::RNGType RNGType;
    //typedef typename Point::FT NT;
    typedef typename Spectrahedron::VT VT;
    unsigned int n =p.dimension();
    RNGType &rng = var.rng;

    Point l = get_direction<RNGType, Point, NT>(n);

    VT pointVT = p.getCoefficients();
    VT lVT = l.getCoefficients();
    std::pair <NT, NT> dbpair = P.boundaryOracle(pointVT, lVT);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
    NT c1 = BoltzmannDirection.dot(b1);
    NT c2 = BoltzmannDirection.dot(b2);

    NT lambda;

    if (c1 > c2) {
        lambda = texp((c1 - c2) / BoltzmannParameter, NT(0), min_plus - max_minus, rng);
        p = b2;
    }
    else {
        lambda = -texp((c2 - c1) / BoltzmannParameter, NT(0), min_plus - max_minus, rng);
        p = b1;
    }

    p = (lambda * l) + p;
}

template <class Point, class Spectrahedron, class Parameters>
void hit_and_run_spec(Point& point,
        Spectrahedron &spectrahedron,
        Parameters &var) {
    typedef typename Spectrahedron::VT VT;
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


//template <class Point, class Parameters>
//void hit_and_run_Boltzmann(Point& point,
//                 Spectrahedron &spectrahedron,
//                 Parameters &var,
//                 const Point& BoltzmannDirection,
//                 double BoltzmannParameter) {
//    typedef typename Parameters::RNGType RNGType;
//    unsigned int n = point.dimension();
//    RNGType &rng = var.rng;
//
//    Point l = get_direction<RNGType, Point, double>(n);
////    VT pointVT = point.getCoefficients();
////    VT lVT = l.getCoefficients();
//    std::pair <double, double> dbpair = spectrahedron.boundaryOracle(point.getCoefficients(), l.getCoefficients());
//    double min_plus = dbpair.first;
//    double max_minus = dbpair.second;
//    Point b1 = (min_plus * l) + point;
//    Point b2 = (max_minus * l) + point;
//    double lambda = boltzmann_distribution(BoltzmannDirection.getCoefficients(), point.getCoefficients(), BoltzmannParameter);
//    point = (lambda * b1);
//    point = ((1 - lambda) * b2) + point;
//}

template <class Point, class Spectrahedron, class Parameters, class VT>
void hit_and_run_spec(Point& point,
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
    std::pair <double, double> dbpair = spectrahedron.boundaryOracleEfficient(pointVT, lVT, a, b);
    double min_plus = dbpair.first;
    double max_minus = dbpair.second;
    Point b1 = (min_plus * l) + point;
    Point b2 = (max_minus * l) + point;
    double lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

template <class Point, class Spectrahedron, class Parameters, class VT, class MT>
void hit_and_run_sampled_covariance_matrix_spec(Point& point,
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
    std::pair <double, double> dbpair = spectrahedron.boundaryOracleEfficient(pointVT, lVT, a, b);
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

template <class Point, class Spectrahedron, class Parameters>
void hit_and_run_spec(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 Point &b1,
                 Point &b2) {
    typedef typename Spectrahedron::VT VT;
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


template <typename NT, class Point, class Spectrahedron, class Parameters, class VT>
void hit_and_run_spec(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 Point &b1,
                 Point &b2,
                 VT& a,
                 NT b) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    std::pair <NT, NT> dbpair = spectrahedron.boundaryOracleEfficient(pointVT, lVT, a, b);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    b1 = (min_plus * l) + point;
    b2 = (max_minus * l) + point;
    NT lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}

template <typename NT, class Point, class Spectrahedron, class SpecSettings, class Parameters>
void hit_and_run_coord_update_spec(Point &p,
                              Spectrahedron &spectrahedron,
                              const unsigned int& rand_coord,
                              const Point& a,
                              const NT & b,
                              Parameters &var,
                              SpecSettings& settings,
                              NT& returnLambda) {
    std::pair<NT, NT> bpair = spectrahedron.boundaryOracleCDHR (p.getCoefficients(), rand_coord, a.getCoefficients(), b, settings);

    if (bpair.first == 0 && bpair.second == 0) {
        returnLambda = 0;
        return;
    }

    typedef typename Parameters::RNGType RNGType;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    RNGType &rng = var.rng;
    NT kapa = urdist(rng);

    returnLambda = bpair.first + kapa * (bpair.second - bpair.first);

    settings.LMIatP.noalias() += returnLambda * spectrahedron.getLMI().getMatrices()[rand_coord];
    p.set_coord(rand_coord, p[rand_coord] + returnLambda);
}

template <class Point, class Spectrahedron, class Parameters, class VT, typename NT, class MT>
void hit_and_run_sampled_covariance_matrix_spec(Point& point,
                 Spectrahedron &spectrahedron,
                 Parameters &var,
                 Point &b1,
                 Point &b2,
                 VT& a,
                 NT b,
                 MT& covarianceMatrix) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = point.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    VT pointVT = point.getCoefficients();
    VT lVT = l.getCoefficients();
    lVT = covarianceMatrix * lVT;
    l = Point(lVT);
    std::pair <NT, NT> dbpair = spectrahedron.boundaryOracleEfficient(pointVT, lVT, a, b);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    b1 = (min_plus * l) + point;
    b2 = (max_minus * l) + point;
    NT lambda = urdist(rng);
    point = (lambda * b1);
    point = ((1 - lambda) * b2) + point;
}


template<class Polytope, class Point, class Parameters, class MT>
void hit_and_run_sampled_covariance_matrix(Point &p,
                                         Polytope &P,
                                         Parameters &var,
                                         Point &b1,
                                         Point &b2,
                                         MT &covarianceMatrix) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    typedef typename Polytope::VT VT;

    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    VT lVT = l.getCoefficients();
    lVT = covarianceMatrix * lVT;
    l = Point(lVT);
    std::pair<NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    b1 = (min_plus * l) + p;
    b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
}

//hit-and-run with random directions and update
// copied this function - need it to return the end points of the segment
template<class Polytope, class Point, class Parameters, class MT>
void hit_and_run_sampled_covariance_matrix(Point &p,
                                           Polytope &P,
                                           Parameters &var,
                                           MT &covarianceMatrix) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    typedef typename Polytope::VT VT;

    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction<RNGType, Point, NT>(n);
    VT lVT = l.getCoefficients();
    lVT = covarianceMatrix * lVT;
    l = Point(lVT);
    std::pair<NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
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
template <class Polytope, class Point, typename NT, class VT>
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

template <class ConvexBody, class Point, class Parameters, typename NT>
void billiard_walk(ConvexBody &P, Point &p, NT che_rad, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev,
                   Parameters &var, bool first = false) {

    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    NT T = urdist(rng) * che_rad;
    Point v = get_direction<RNGType, Point, NT>(n);
    Point p0=p, v_prev=v;
    std::pair<NT, int> pbpair;
    int it=0;

    NT lambda_max=0.0;
    if (first) {

        //p.print();
        //std::cout<<"[first][in bill] is_in = "<<P.is_in(p)<<std::endl;
        //std::cout<<"[1]before line intersection"<<std::endl;
        pbpair = P.line_positive_intersect(p, v, Ar, Av);
        //while(pbpair.first < 0.00001 && pbpair.second != P.num_of_hyperplanes()){
        // /   std::cout<<"lambda_prev = "<<pbpair.first<<std::endl;
        //    v = get_direction<RNGType, Point, NT>(n);
        //    //p = ((-0.01*lambda_prev)*v)+p;
        //    pbpair = P.line_positive_intersect(p, v, Ar, Av);
        //    std::cout<<"lambda_meta = "<<pbpair.first<<std::endl;
        //    sleep(3);
        //}
        //std::cout<<"[1]after line intersection"<<std::endl;
        if (T <= pbpair.first) {
            p = (T * v) + p;
            lambda_prev = T;
            return;
        }
        lambda_prev = 0.995 * pbpair.first;
        //std::cout<<"[first]lambda_pos ="<<lambda_prev<<" T = "<<T<<std::endl;
        p = (lambda_prev * v) + p;
        T -= lambda_prev;
        //std::cout<<"[1]before reflection"<<std::endl;
        v_prev=v;
        P.compute_reflection(v, p, pbpair.second);
        //if(P.is_in(p)==0) throw "Point out!";
        //if (lambda_prev<0.0000000001) throw "small lambda!";
        //std::cout<<"[1]after reflection"<<std::endl;
    }

    while (it<3*n) {

        //p.print();
        //std::cout<<"bef [in bill] is_in = "<<P.is_in(p)<<std::endl;
        //std::cout<<"[2]before line intersection"<<std::endl;
        pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev);
        //while(pbpair.first < 0.00001 && pbpair.second != P.num_of_hyperplanes()){
        //std::cout<<"lambda_prev = "<<pbpair.first<<std::endl;
        //p = ((-0.01*lambda_prev)*v_prev)+p;
        //pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev);
        //std::cout<<"lambda_meta = "<<pbpair.first<<std::endl;
        //sleep(3);
        //}

        //std::cout<<"[2]after line intersection"<<std::endl;
        //std::cout<<"lambda_pos ="<<pbpair.first<<" T = "<<T<<std::endl;
        if (T <= pbpair.first) {
            p = (T * v) + p;
            lambda_prev = T;
            break;
        }

        lambda_prev = 0.995 * pbpair.first;
        //std::cout<<"lambda_pos ="<<lambda_prev<<" T = "<<T<<std::endl;
        p = (lambda_prev * v) + p;

        T -= lambda_prev;
        //std::cout<<"[2]before reflection"<<std::endl;
        v_prev=v;
        P.compute_reflection(v, p, pbpair.second);
        //std::cout<<"[in bill] is_in = "<<P.is_in(p)<<std::endl;
        //if(P.is_in(p)==0) throw "Point out!";
        //if (lambda_prev<0.00000000001) throw "small lambda!";
        it++;
        //std::cout<<"[2]after reflection"<<std::endl;
    }
//    if(it==10*n){
        //for (int i = 0; i < 1; ++i) std::cout<<"IT IS 10*N!!!!!!!!"<<std::endl;
        //for (int i = 0; i < 1; ++i) std::cout<<"lambda_pos ="<<lambda_prev<<" T = "<<T<<std::endl;
        //std::cout<<"first = "<<first<<std::endl;
        //sleep(3);
        //throw "Point out!";
//        p=p0;
//    }
}


//template <class Point, class Parameters, typename NT>
//void billiard_walk(Spectrahedron &spectrahedron, Point &p, NT che_rad, MT& LMIatP, Parameters &var, const Point& a, double b, std::list<Point>& randPoints, bool first = true) {
//
//    typedef typename Parameters::RNGType RNGType;
//    unsigned int n = spectrahedron.getLMI().getDim();
//    RNGType &rng = var.rng;
//    boost::random::uniform_real_distribution<> urdist(0, 1);
//    NT T = urdist(rng) * che_rad;
//    Point v = get_direction<RNGType, Point, NT>(n);
//    int it = 0;
//    std::pair<double, bool> pair;
//    NT lambda_max = 0.0;
//    MT B;
//    VT genEivector;
////    VT pVT = p.getCoefficients();
//
////    std::pair<double, double> bpair;
//
//    if (first) {
//
//        pair = spectrahedron.boundaryOraclePositive(p.getCoefficients(), v.getCoefficients(), a.getCoefficients(), b, LMIatP, B, genEivector, true);
////        bpair = spectrahedron.boundaryOracleEfficient(pVT, v, a, b);
//        double lambda = pair.first;
//
////        std::cout<<pair.first << " lambda " <<std::endl;
//        lambda = 0.995 * lambda;
//        if (lambda == 0)
//            return;
//
//        if (T <= lambda) {
//            p = (T * v) + p;
//            LMIatP += T*B;
////            std::cout<< " =========\n\n\n\n\n\n" <<std::endl;
//            return;
//        }
//
//        p = (lambda * v) + p;
//        randPoints.push_back(p);
//        LMIatP += lambda*B;
//        T -= lambda;
//
//        if (pair.second) {
//            //we hit the cutting plane
//            Point reflection = ((-2.0 * v.dot(a)) * a);
//            v = v + reflection;
//        }
//        else
//            spectrahedron.compute_reflection(genEivector, v);
//    }
//
//    while (it < 10 * n) {
//
//        pair = spectrahedron.boundaryOraclePositive(p.getCoefficients(), v.getCoefficients(), a.getCoefficients(), b, LMIatP, B, genEivector, false);
//
//        double lambda = pair.first;
////        std::cout<<pair.first << " lambda " <<std::endl;
//        lambda = 0.995 * lambda;
//        if (lambda == 0) {
//            return;
//        }
//
//        if (T <= lambda) {
//            p = (T * v) + p;
//            LMIatP += T*B;
//            break;
//        }
//
//        p = (lambda * v) + p;
//        randPoints.push_back(p);
//
//        LMIatP += lambda*B;
//        T -= lambda;
//
//        if (pair.second) {
//            //we hit the cutting plane
//            Point reflection = ((-2.0 * v.dot(a)) * a);
//            v = v + reflection;
//        }
//        else
//            spectrahedron.compute_reflection(genEivector, v);
//
//        it++;
//    }
//}

//template <class Point, class Parameters, typename NT>
//void billiard_walk(Spectrahedron &spectrahedron, Point &p, NT che_rad, MT& LMIatP, Parameters &var, const Point& a, double b, Point& v, bool first = true) {
//
//    typedef typename Parameters::RNGType RNGType;
//    unsigned int n = spectrahedron.getLMI().getDim();
//    RNGType &rng = var.rng;
//    boost::random::uniform_real_distribution<> urdist(0, 1);
//    NT T = urdist(rng) * che_rad;
//    int it = 0;
//    std::pair<double, bool> pair;
//    NT lambda_max = 0.0;
//    MT B;
//    VT genEivector;
////    VT pVT = p.getCoefficients();
//
////    std::pair<double, double> bpair;
//
//    if (first) {
//
//        pair = spectrahedron.boundaryOraclePositive(p.getCoefficients(), v.getCoefficients(), a.getCoefficients(), b, LMIatP, B, genEivector, true);
////        bpair = spectrahedron.boundaryOracleEfficient(pVT, v, a, b);
//        double lambda = pair.first;
//
////        std::cout<<pair.first << " vs " << bpair.first<<std::endl;
//        lambda = 0.995 * lambda;
//        if (lambda == 0)
//            return;
//
//        if (T <= lambda) {
//            p = (T * v) + p;
//            LMIatP += T*B;
////            std::cout<< " =========\n\n\n\n\n\n" <<std::endl;
//            return;
//        }
//
//        p = (lambda * v) + p;
//        LMIatP += lambda*B;
//        T -= lambda;
//
//        if (pair.second) {
//            //we hit the cutting plane
//            Point reflection = ((-2.0 * v.dot(a)) * a);
//            v = v + reflection;
//        }
//        else
//            spectrahedron.compute_reflection(genEivector, v);
//    }
//
//    while (it < 10 * n) {
//
//        pair = spectrahedron.boundaryOraclePositive(p.getCoefficients(), v.getCoefficients(), a.getCoefficients(), b, LMIatP, B, genEivector, false);
//
//        double lambda = pair.first;
//
//        lambda = 0.995 * lambda;
//        if (lambda == 0) {
//            return;
//        }
//
//        if (T <= lambda) {
//            p = (T * v) + p;
//            LMIatP += T*B;
//            break;
//        }
//
//        p = (lambda * v) + p;
//
//        LMIatP += lambda*B;
//        T -= lambda;
//
//        if (pair.second) {
//            //we hit the cutting plane
//            Point reflection = ((-2.0 * v.dot(a)) * a);
//            v = v + reflection;
//        }
//        else
//            spectrahedron.compute_reflection(genEivector, v);
//
//        it++;
//    }
//}


template <class Spectrahedron, class Point, class Parameters, class SpecSettings, typename NT>
double billiard_walk(Spectrahedron &spectrahedron, Point &p, const NT& che_rad, Parameters &var, const Point& a, const NT& b,
                         SpecSettings& settings, bool shake_and_bake_directions = false) {

    typedef typename Parameters::RNGType RNGType;
    unsigned int n = spectrahedron.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    NT T = urdist(rng) * che_rad;
    Point v, p0 = p;

    if (!shake_and_bake_directions)
        v = get_direction<RNGType, Point, NT>(n);
    else
        rsab_get_direction<RNGType, Point, NT, Parameters> (n,a,v,var);

    NT it = 0;
    std::pair<NT, bool> pair;
    NT factor = 100.0;

    while (it < factor*n) {

        pair = spectrahedron.boundaryOracleBilliard(p.getCoefficients(), v.getCoefficients(), a.getCoefficients(), b, settings, it != 0 && shake_and_bake_directions);

        NT lambda = pair.first;

        if (lambda == 0) {
            return it;
        }

        it++;

        lambda = 0.995 * lambda;

        if (T <= lambda) {
            p += (T * v);
            settings.LMIatP.noalias() += T*settings.B;
            break;
        }

        p += (lambda * v);
        settings.LMIatP.noalias() += lambda*settings.B;
        T -= lambda;

        if (it >= factor*n)
            break;

        if (pair.second) {
            //we hit the cutting plane
            double l = -2.0 * v.dot(a);
            v +=  l * a;
            settings.B += l * settings.Obj;
            settings.computeB = false;
        }
        else
            spectrahedron.compute_reflection(settings.genEigenvector, v, p);

    }

    if (it>=factor*n) {
        std::cout<<"limit reached"<<std::endl;
        p=p0;
    }

    return it;
}

template <class Spectrahedron, class Point, class Parameters, class SpecSettings, typename NT>
void HMC_boltzmann_reflections(Spectrahedron &spectrahedron, Point &p, NT che_rad, Parameters &var,  Point& objectiveFunction,  NT& temperature,
        SpecSettings& settings) {

    typedef typename Parameters::RNGType RNGType;
    typedef typename Spectrahedron::MT MT;
    typedef typename Spectrahedron::VT VT;
    unsigned int n = spectrahedron.getLMI().getDim();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    NT T = urdist(rng) * che_rad;//0.7*che_rad;
    Point v(n);
    int it = 0, max_it = n;
    NT lambda_max = 0.0;
    VT pVT = p.get_coefficients();

    if(n==2 || n==3) max_it=10*n;

//    Point p1 = p;

//    temperature = 10;
//    while (temperature > 0.0001) {

//        pVT(1) = 0;
//        pVT(0) = -0.3;
//        T = urdist(rng) * che_rad;
//        it = 0;

//        std::cout << "TEMP " << temperature << "\n";
//        v(0) = 1;
//        v(1) = 0.5;

        //if (settings.isBoundaryPoint && !settings.first) {
      //      rsab_get_direction<RNGType, Point, NT, Parameters>(n, settings.gradient, v, var);
        //}
        //else
        v = get_direction<RNGType, Point, NT>(n);

//        if (first) {
//            std::cout << pVT.transpose() << "\n";
//            std::cout << v.transpose() << "\n";
//            double lambda = spectrahedron.boundaryOracle_Boltzmann_HMC(pVT, v, objectiveFunction, temperature, LMIatP,
//                                                                       B1, B2, genEivector, true);


//            double lambdaSafe = 0.99 * lambda;

//            MT C = LMIatP + (lambdaSafe * B1) + (((-0.5 * lambdaSafe * lambdaSafe) / temperature) * B2);


//            if (lambda == 0) {
//                return;
//                temperature = temperature / 2;
//                continue;
//            }

//            if (T <= lambdaSafe) {
//                pVT = (-0.5 / temperature) * objectiveFunction * T * T + v * T + pVT;
//                p = Point(pVT);
//                Point ss(objectiveFunction);
//                std::cout << " $ending " << p.dot(ss) << "\n";
//                LMIatP += (T * B1) + (((-0.5 * T * T) / temperature) * B2);
//                return;
//                temperature = temperature /2;
//                continue;
//            }

//            pVT = (-0.5 / temperature) * objectiveFunction * lambdaSafe * lambdaSafe + v * lambdaSafe + pVT;
//            LMIatP += (lambdaSafe * B1) + (((-0.5 * lambdaSafe * lambdaSafe) / temperature) * B2);
//            T -= lambdaSafe;

//            v = -objectiveFunction * (lambda / (T)) + v;
//            std::cout << v.transpose() << "\n";
//            spectrahedron.compute_reflection(genEivector, v, C);
//            std::cout << v.transpose() << "\n";
//        }


        while (it < max_it) {
//            Point oo(pVT);
            typedef std::numeric_limits<double> dbl;

            std::cout.precision(dbl::max_digits10);
//            std::cout << " $ " << oo.dot(ss) << " " << temperature << " " << che_rad << "\n";
//            std::cout << pVT.transpose() << "\n";
//            v.print();

//            p.print();
//            v.print();
//
//            if (!settings.first)
//            settings.gradient.print();
            NT lambda = spectrahedron.boundaryOracle_Boltzmann_HMC(p, v, objectiveFunction, temperature, settings);
            settings.first = false;

            MT C;// = LMIatP + (lambda * B1) + (((-0.5 * lambda * lambda) / temperature) * B2);

            NT lambdaSafe = 0.99 * lambda;

            settings.isBoundaryPoint = lambdaSafe < settings.epsilon;

            if (lambda == 0) {
//                p = Point(pVT);
                return;
            }

            if (T <= lambdaSafe) {
                Point t = (-0.5 / temperature) * objectiveFunction * T * T;
                p = p + t;
                t = v * T;
                p = p + t;
                settings.B0 += (T * settings.B1) + (((-0.5 * T * T) / temperature) * settings.B2);
                break;
            }

            Point t = (-0.5 / temperature) * objectiveFunction * lambdaSafe * lambdaSafe;
            p = p + t;
            t = v * lambdaSafe;
            p = p + t;

            settings.B0 += (lambda * settings.B1) + (((-0.5 * lambdaSafe * lambdaSafe) / temperature) * settings.B2);
            T -= lambdaSafe;

            v = -1*objectiveFunction * (lambda / (temperature)) + v;
//            std::cout << v.transpose() << "\n";
            spectrahedron.compute_reflection(settings, v);

            it++;
        }
//        temperature = temperature /2;
//    }
//    exit(0);
//    std::cout <<"end\n";
      //p = Point(pVT);
}



#endif //RANDOM_SAMPLERS_H
