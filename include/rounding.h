// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// The functions in this header file call Bojan Nikolic <bojan@bnikolic.co.uk> implementation
// of Todd and Yildirim algorithm in "On Khachiyan's Algorithm for the Computation of Minimum Volume Enclosing Ellipsoids", 2005
// Initial version 2010, part of BNMin1 and is licensed under GNU General Public License version 2

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ROUNDING_H
#define ROUNDING_H


#include "khach.h"


// ----- ROUNDING ------ //
// main rounding function
template <typename Polytope, typename Point, typename Parameters, typename NT>
std::pair <NT, NT> rounding_min_ellipsoid(Polytope &P , const std::pair<Point,NT> &InnerBall, const Parameters &var) {

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Parameters::RNGType RNGType;

    unsigned int n=var.n, walk_len=var.walk_steps, i, j = 0;
    Point c = InnerBall.first;
    NT radius = InnerBall.second;
    std::list<Point> randPoints; //ds for storing rand points
    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        // If P is not a V-Polytope or number_of_vertices>20*domension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        p = p + c;

        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 10*n, randPoints, var);
        // 3. Sample points from P
        unsigned int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
        randPoints.clear();
        if (var.bill_walk) {
            rand_point_generator(P, p, num_of_samples, 5, randPoints, var);
        } else {
            rand_point_generator(P, p, num_of_samples, 10 + n / 10, randPoints, var);
        }
    }

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
    typename std::list<Point>::iterator rpit=randPoints.begin();

    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        for (i=0 ; i<rpit->dimension(); i++){
            Ap(i,j)=double((*rpit)[i]);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n); //TODO: remove dependence on ublas and copy to eigen
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

    MT E(n,n);
    VT e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(unsigned int i=0; i<n; i++){
        e(i)=NT(c2(i));
        for (unsigned int j=0; j<n; j++){
            E(i,j)=NT(Q(i,j));
        }
    }


    //Find the smallest and the largest axes of the elliposoid
    Eigen::EigenSolver<MT> eigensolver(E);
    NT rel = std::real(eigensolver.eigenvalues()[0]);
    NT Rel = std::real(eigensolver.eigenvalues()[0]);
    for(unsigned int i=1; i<n; i++){
        if(std::real(eigensolver.eigenvalues()[i])<rel) rel=std::real(eigensolver.eigenvalues()[i]);
        if(std::real(eigensolver.eigenvalues()[i])>Rel) Rel=std::real(eigensolver.eigenvalues()[i]);
    }

    Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
    MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    P.shift(e);

    MT L_1 = L.inverse();
    P.linear_transformIt(L_1.transpose());

    return std::pair<NT, NT> (L_1.determinant(),rel/Rel);
}


// ----- ROUNDING ------ //
// main rounding function
template <typename MT, typename VT, typename Polytope, typename Point, typename Parameters, typename NT>
std::pair< std::pair<MT, VT>, NT >  rounding_min_ellipsoid(Polytope &P , const std::pair<Point,NT> &InnerBall,
                                         const Parameters &var) {

    //typedef typename Polytope::VT 	VT;
    typedef typename Parameters::RNGType RNGType;

    unsigned int n=var.n, walk_len=var.walk_steps, i, j = 0;
    Point c = InnerBall.first;
    NT radius = InnerBall.second;
    std::list<Point> randPoints; //ds for storing rand points
    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        // If P is not a V-Polytope or number_of_vertices>20*domension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        p = p + c;

        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 10*n, randPoints, var);
        // 3. Sample points from P
        unsigned int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
        randPoints.clear();
        if (var.bill_walk) {
            rand_point_generator(P, p, num_of_samples, 5, randPoints, var);
        } else {
            rand_point_generator(P, p, num_of_samples, 10 + n * 5, randPoints, var);
        }
    }

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
    typename std::list<Point>::iterator rpit=randPoints.begin();

    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        for (i=0 ; i<rpit->dimension(); i++){
            Ap(i,j)=double((*rpit)[i]);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n); //TODO: remove dependence on ublas and copy to eigen
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

    MT E(n,n);
    VT e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(unsigned int i=0; i<n; i++){
        e(i)=NT(c2(i));
        for (unsigned int j=0; j<n; j++){
            E(i,j)=NT(Q(i,j));
        }
    }


    //Find the smallest and the largest axes of the elliposoid
    Eigen::EigenSolver<MT> eigensolver(E);
    NT rel = std::real(eigensolver.eigenvalues()[0]);
    NT Rel = std::real(eigensolver.eigenvalues()[0]);
    for(unsigned int i=1; i<n; i++){
        if(std::real(eigensolver.eigenvalues()[i])<rel) rel=std::real(eigensolver.eigenvalues()[i]);
        if(std::real(eigensolver.eigenvalues()[i])>Rel) Rel=std::real(eigensolver.eigenvalues()[i]);
    }

    Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
    MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    P.shift(e);

    MT L_1 = L.inverse();
    P.linear_transformIt(L_1.transpose());

    return std::pair< std::pair<MT, VT>, NT > (std::pair<MT, VT>(L_1, e), L_1.determinant());
}


template <typename Polytope>
void get_vpoly_center(Polytope &P) {

    typedef typename Polytope::NT 	NT;
    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Polytope::PolytopePoint 	Point;

    unsigned int n = P.dimension();

    std::list<Point> randPoints; //ds for storing rand points
    if (!P.get_points_for_rounding(randPoints)) {
        P.shift(P.get_mean_of_vertices().getCoefficients());
    } else {

        boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
        typename std::list<Point>::iterator rpit=randPoints.begin();
        for (int j=0 ; rpit!=randPoints.end(); rpit++, j++) {
            for (int i=0 ; i<rpit->dimension(); i++){
                Ap(i,j)=double((*rpit)[i]);
            }
        }
        boost::numeric::ublas::matrix<double> Q(n,n);
        boost::numeric::ublas::vector<double> c2(n);
        size_t w=1000;
        KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

        VT e(n);
        //Get ellipsoid matrix and center as Eigen objects
        for(unsigned int i=0; i<n; i++){
            e(i)=NT(c2(i));
        }
        P.shift(e);
    }

}


#endif
