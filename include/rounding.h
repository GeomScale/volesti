// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ROUNDING_H
#define ROUNDING_H


#include "khach.h"

/* EXPERIMENTAL
//ROUNDING
template <class MT, class Point>
MT getPointsMat(std::list<Point> randPoints, int dim){

    MT S(randPoints.size(),dim);
    for(unsigned int i=0; i<randPoints.size(); i++){
        Point p=randPoints.front();
        randPoints.pop_front();
        for (unsigned int j=0; j<dim; j++){
            S(i,j)=p[j];
        }
    }
    
    return S;
}


template <typename NT, class Polytope, class Point, class Parameters>
std::pair<Point, NT> approx_R(Polytope &P, Parameters var){
    std::pair<Point,NT> InnerBall=P.ComputeInnerBall();
    Point c=InnerBall.first;
    NT radius = InnerBall.second;

    int n=var.n, walk_len=var.walk_steps;
    //Random_points_on_sphere_d<Point> gen (n, radius);
    //Point p = gen.sample_point(var.rng);
    Point p = get_point_on_Dsphere(n, radius);
    p = p + c;
    std::list<Point> randPoints; //ds for storing rand points
    //use a large walk length e.g. 1000
    rand_point_generator(P, p, 1, 50*n, randPoints, var);
    //if (print) std::cout<<"First random point: "<<p<<std::endl;

    // 3. Sample points from P
    //randPoints.push_front(p);
    int num_of_samples = std::pow(1.0,-2) * 400 * n * std::log(n);;//this is the number of sample points will used to compute min_ellipoid
    //if(print) std::cout<<"\nCompute "<<num_of_samples<<" random points in P"<<std::endl;
    rand_point_generator(P, p, num_of_samples*10, 1, randPoints, var);
    NT current_dist, max_dist;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    return std::pair<Point,NT> (c,R);
}


//Needs developing (experimental)
template <class Polytope, class Point, class Parameters, typename  NT>
NT rounding_SVD(Polytope &P , Point c, NT radius, Parameters &var){
    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    int n=var.n, walk_len=var.walk_steps;
    bool print=var.verbose;
    // 1. Compute the Chebychev ball (largest inscribed ball) with center and radius 
	//Point c(n);       //center
    //K radius;
    //P.chebyshev_center(c,radius);
    
    // 2. Generate the first random point in P
  // Perform random walk on random point in the Chebychev ball 
	//if (print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
	//vars var2=var;
	//var2.coordinate=false;
    Point p = get_point_on_Dsphere(n, radius);
	p = p + c;
	std::list<Point> randPoints; //ds for storing rand points
	//use a large walk length e.g. 1000
	rand_point_generator(P, p, 1, 50*n, randPoints, var);
	//if (print) std::cout<<"First random point: "<<p<<std::endl;
    // 3. Sample points from P
	//randPoints.push_front(p);
	int num_of_samples = std::pow(1.0,-2) * 400 * n * std::log(n);;//this is the number of sample points will used to compute min_ellipoid
	//if(print) std::cout<<"\nCompute "<<num_of_samples<<" random points in P"<<std::endl;
	rand_point_generator(P, p, num_of_samples*10, 1, randPoints, var);
    NT current_dist, max_dist;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    if(print) std::cout<<"R = "<<max_dist<<" r = "<<radius<<"ratio R/r = "<<R<<"\n"<<std::endl;
    
    // 4. Compute the transformation matrix T
    MT T = Eigen::MatrixXd::Identity(n,n);
    bool well_rounded=false;
    int t=8*n*n*n;
    //int t=var.m;
    int tries=0;
    MT S=Eigen::MatrixXd::Identity(n,n);
    std::pair<Point,NT> res;
    while(!well_rounded){
        tries++;
        randPoints.clear();
        Polytope P2(P);
        P2.linear_transformIt(T.inverse());   //We have to sample from the transformed body
        res=solveLP(P2.get_matrix(), P2.dimension());
        c=res.first;
        p = get_point_on_Dsphere(n, res.second);
        p = p + c;
        rand_point_generator(P2, p, 1, 50*n, randPoints, var);
        randPoints.clear();
        rand_point_generator(P2, p, t, 1, randPoints, var);
        MT PM=getPointsMat<MT>(randPoints,n);
        Eigen::JacobiSVD<MT> svd(PM, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if(print) std::cout<<svd.singularValues()<<"\n"<<std::endl;
        NT min=svd.singularValues()(0);
        for(unsigned int i=1; i<n; i++){
            if(svd.singularValues()(i)<min){
                min=svd.singularValues()(i);
            }
        }
        for(unsigned int i=0; i<n; i++){
            S(i,i)=svd.singularValues()(i)/min;
        }
        if(print) std::cout<<S<<"\n"<<std::endl;

        T=svd.matrixV()*S.inverse()*T;
        if(print) std::cout<<T<<"\n"<<std::endl;
        well_rounded=true;
        for (unsigned int i=0; i<n; i++){
            if (S(i,i)>2.0){
                if (tries>((int)log2(R))){
                    t=t*2;
                    tries=0;
                }
                well_rounded=false;
                break;
            }
        }
        if (well_rounded){
            P.linear_transformIt(T.inverse());
        }
    }
    
    return std::abs(T.determinant());
}*/


// ----- ROUNDING ------ //
// main rounding function
template <class Polytope, class Point, class Parameters, typename NT>
std::pair <NT, NT> rounding_min_ellipsoid(Polytope &P , std::pair<Point,NT> InnerBall, Parameters &var) {

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
    typename std::vector<NT>::iterator qit;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            Ap(i,j)=double(*qit);
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


// -------- ROTATION ---------- //
/*
template <class T>
double rotating_old(T &P){

  bool print = true; 
  if(print) std::cout<<"\nRotate..."<<std::endl;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MT;
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> VT;
  
  int m=P.num_of_hyperplanes();
  int d=P.dimension();
  
  MT A(m,d);
  VT b(m);
  for(int i=0; i<m; ++i){
		b(i) = P.get_coeff(i,0);
		for(int j=1; j<d+1; ++j){
		  A(i,j-1) = P.get_coeff(i,j);
		}
	}
  std::cout<<A<<"\n"<<b<<std::endl;
  
  // Construct rotation matrix by 30o
  MT R(m,m);
  for(int i=0; i<m; ++i){
		for(int j=0; j<m; ++j){
		  R(i,j) = 0;
		  if (i==j) R(i,j) = 1;
		  if (i==0 && j==0) R(i,j) = std::sqrt(3)/2;
		  if (i==1 && j==0) R(i,j) = -1.0/2.0;
		  if (i==0 && j==1) R(i,j) = 1.0/2.0;
		  if (i==1 && j==1) R(i,j) = std::sqrt(3)/2;
		  if (i==d && j==d) R(i,j) = std::sqrt(3)/2;
		  if (i==d+1 && j==d) R(i,j) = -1.0/2.0;
		  if (i==d && j==d+1) R(i,j) = 1.0/2.0;
		  if (i==d+1 && j==d+1) R(i,j) = std::sqrt(3)/2;
		}
	}
	//std::cout<<R<<std::endl;
	
	A = R*A;
	//b = R*b;
	
	std::cout<<A<<"\n"<<b<<std::endl;
	
	// Write changes (actually perform rotation) to the polytope!
	for(int i=0; i<m; ++i){
		P.put_coeff(i,0,b(i));
		for(int j=1; j<d+1; ++j){
		  P.put_coeff(i,j,A(i,j-1));
		}
	}
	
  std::cout<<R.determinant()<<"\n"<<b<<std::endl;
  
	return R.determinant();
}*/

// -------- ROTATION ---------- //
template <class MT, class Polytope>
MT rotating(Polytope &P){

    typedef boost::mt19937    RNGType;
    //typedef typename Polytope::MT 	MT;

    boost::random::uniform_real_distribution<> urdist(-1.0, 1.0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    unsigned int n = P.dimension();

    // pick a random rotation
    MT R(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            R(i,j) = urdist(rng);
        }
    }

    Eigen::JacobiSVD<MT> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // apply rotation to the polytope P
    P.linear_transformIt(svd.matrixU());

    return svd.matrixU().inverse();
}

#endif
