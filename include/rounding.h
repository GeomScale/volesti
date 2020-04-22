// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ROUNDING_H
#define ROUNDING_H


#include "khach.h"
//#include "Eigen"

// ----- ROUNDING ------ //
// main rounding function
template <class Polytope, class Point, class Parameters, typename NT>
std::pair <NT, NT> rounding_min_ellipsoid(Polytope &P , std::pair<Point,NT> InnerBall, Parameters &var) {

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Parameters::RNGType RNGType;
    typedef Eigen::Matrix<NT, Eigen::Dynamic,1> Coeff;

    unsigned int n=var.n, walk_len=var.walk_steps, i, j = 0;
    Point c = InnerBall.first;
    NT radius = InnerBall.second;
    std::list<Point> randPoints; //ds for storing rand points
    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        // If P is not a V-Polytope or number_of_vertices>20*domension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_on_Dsphere<RNGType, Point>(n, radius);
        p = p + c;

        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 50*n, randPoints, var);
        // 3. Sample points from P
        unsigned int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
        randPoints.clear();
        rand_point_generator(P, p, num_of_samples, walk_len, randPoints, var);
        /*NT current_dist, max_dist;
        for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
            current_dist=(*pit-c).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
        R=max_dist/radius;*/
    }

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
    typename std::list<Point>::iterator rpit=randPoints.begin();


    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        Coeff coeffs = rpit->getCoefficients();
        for (i=0 ; i<rpit->dimension(); i++){
            Ap(i,j)=double(coeffs(i));
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n);
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

template <class Spectrahedron, class Point, class Parameters, class SpecSettings, typename NT>
void preproccess_spectrahedron(Spectrahedron &SP, Point &p, Parameters &var, SpecSettings settings,
                                   NT &rand_value, NT &diam, NT &radius, bool &rounding) {

    //typedef typename Polytope::NT 	NT;
    typedef typename Spectrahedron::MT MT;
    typedef typename Spectrahedron::VT VT;
    //typedef typename Polytope::PolytopePoint 	Point;

    unsigned int n = var.n, walk_len = 1, i, j;
    std::list <Point> randPoints;
    rand_value = 1.0;
    //var.cdhr_walk = false;
    MT E(n, n);
    VT e(n);
    boost::numeric::ublas::matrix<double> Ap(n, 10 * n);
    NT max_diam = 0.0, diam_iter, ratio1 = 0.0;

    //std::cout<<"computing enclosed ball.."<<std::endl;
    SP.ComputeInnerBall(diam, radius);
    //std::cout<<"enclosed ball computed.."<<std::endl;
    var.diameter = diam;
    var.che_rad = radius;

    int count = 0;

    while(true) {

        count++;

        randPoints.clear();
        //std::cout<<"Sampling 10d points from P.."<<std::endl;
        rand_point_generator_spec(SP, p, 10 * n, 1, randPoints, var, settings);
        //std::cout<<"points sampled.."<<std::endl;
        //boundary_rand_point_generator(P, p, 50*n, walk_len, randPoints, var);

        typename std::list<Point>::iterator rpit = randPoints.begin();
        //typename std::vector<NT>::iterator qit;
        VT qq(n);
        j = 0, i = 0;
        for (; rpit != randPoints.end(); rpit++, j++) {
            qq = (*rpit).get_coefficients();
            //i = 0;
            for (i=0; i<n; i++) {
                Ap(i, j) = double(qq(i));
            }
        }
        boost::numeric::ublas::matrix<double> Q(n, n);
        boost::numeric::ublas::vector<double> c2(n);
        size_t w = 1000;
        KhachiyanAlgo(Ap, 0.01, w, Q, c2); // call Khachiyan algorithm



        //Get ellipsoid matrix and center as Eigen objects
        for (unsigned int i = 0; i < n; i++) {
            e(i) = NT(c2(i));
            for (unsigned int j = 0; j < n; j++) {
                E(i, j) = NT(Q(i, j));
            }
        }


        if(!rounding && count==1) break;

        if (rounding || count<=1) {
            Eigen::LLT <MT> lltOfA(E); // compute the Cholesky decomposition of E
            MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
            MT L_1 = L.inverse();

            rand_value *= L_1.determinant();
            SP.shift(e);
            SP.linear_transformIt(L_1.transpose());
            SP.ComputeInnerBall(diam, radius);
            var.che_rad = radius;
            var.diameter = diam;
            p = Point(n);
            SP.set_LMIatP_A0(settings);
            rounding = false;
        } else {
            break;
        }

    }
    SP.shift(e);
    settings.LMIatP = SP.getLMI().getA0();

    SP.ComputeInnerBall(diam, radius);
    var.che_rad = radius;
    var.diameter = diam;
    //std::cout<<"preproccess completed.."<<std::endl;

}

#endif
