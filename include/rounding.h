// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ROUNDING_H
#define ROUNDING_H

//ROUNDING

Eigen::MatrixXd getPointsMat(std::list<Point> randPoints, int dim){
    Eigen::MatrixXd S(randPoints.size(),dim);
    for(int i=0; i<randPoints.size(); i++){
        Point p=randPoints.front();
        randPoints.pop_front();
        for (int j=0; j<dim; j++){
            S(i,j)=p[j];
        }
    }
    
    return S;
}


template <class T1>
std::pair<Point, NT> approx_R(T1 &P, vars var){
    std::pair<Point,double> Cheb_ball=solveLP(P.get_matrix(), P.dimension());
    Point c=Cheb_ball.first;
    NT radius = Cheb_ball.second;

    int n=var.n, walk_len=var.walk_steps;
    Random_points_on_sphere_d<Point> gen (n, radius);
    Point p = gen.sample_point(var.rng);
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
    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    return std::pair<Point,NT> (c,R);
}


//Needs developing
template <class T1>
NT rounding_SVD(T1 &P , Point c, NT radius, vars &var){
    //typedef typename T1::FT 	K;
    int n=var.n, walk_len=var.walk_steps;
    bool print=var.verbose;
    // 1. Compute the Chebychev ball (largest inscribed ball) with center and radius 
	//Point c(n);       //center
    //K radius;
    //P.chebyshev_center(c,radius);
    
    // 2. Generate the first random point in P
  // Perform random walk on random point in the Chebychev ball 
	//if (print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
	vars var2=var;
	var2.coordinate=false;
	Random_points_on_sphere_d<Point> gen (n, radius);
	Point p = gen.sample_point(var.rng);
	p = p + c;
	std::list<Point> randPoints; //ds for storing rand points
	//use a large walk length e.g. 1000
	rand_point_generator(P, p, 1, 50*n, randPoints, var2);
	//if (print) std::cout<<"First random point: "<<p<<std::endl;
    // 3. Sample points from P
	//randPoints.push_front(p);
	int num_of_samples = std::pow(1.0,-2) * 400 * n * std::log(n);;//this is the number of sample points will used to compute min_ellipoid
	//if(print) std::cout<<"\nCompute "<<num_of_samples<<" random points in P"<<std::endl;
	rand_point_generator(P, p, num_of_samples*10, 1, randPoints, var2);
    NT current_dist, max_dist;
    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    if(print) std::cout<<"R = "<<max_dist<<" r = "<<radius<<"ratio R/r = "<<R<<"\n"<<std::endl;
    
    // 4. Compute the transformation matrix T
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(n,n);
    bool well_rounded=false;
    int t=8*n*n*n;
    //int t=var.m;
    int tries=0;
    Eigen::MatrixXd S=Eigen::MatrixXd::Identity(n,n);
    std::pair<Point,NT> res;
    while(!well_rounded){
        tries++;
        randPoints.clear();
        T1 P2(P);
        P2.linear_transformIt(T.inverse());   //We have to sample from the transformed body
        res=solveLP(P2.get_matrix(), P2.dimension());
        c=res.first;
        Random_points_on_sphere_d<Point> gen (n, res.second);
        p = gen.sample_point(var.rng);
        p = p + c;
        rand_point_generator(P2, p, 1, 50*n, randPoints, var2);
        randPoints.clear();
        rand_point_generator(P2, p, t, 1, randPoints, var2);
        Eigen::MatrixXd PM=getPointsMat(randPoints,n);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(PM, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if(print) std::cout<<svd.singularValues()<<"\n"<<std::endl;
        NT min=svd.singularValues()(0);
        for(int i=1; i<n; i++){
            if(svd.singularValues()(i)<min){
                min=svd.singularValues()(i);
            }
        }
        for(int i=0; i<n; i++){
            S(i,i)=svd.singularValues()(i)/min;
        }
        if(print) std::cout<<S<<"\n"<<std::endl;

        T=svd.matrixV()*S.inverse()*T;
        if(print) std::cout<<T<<"\n"<<std::endl;
        well_rounded=true;
        for (int i=0; i<n; i++){
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
}


// ----- ROUNDING ------ //
template <class T1>
std::pair<double,double> rounding_min_ellipsoid(T1 &P , std::pair<Point,double> CheBall, vars &var){
    int n=var.n, walk_len=var.walk_steps;
    bool print=var.verbose;
    Point c = CheBall.first;
    NT radius = CheBall.second;
    // 2. Generate the first random point in P
    // Perform random walk on random point in the Chebychev ball
    Random_points_on_sphere_d<Point> gen (n, radius);
    Point p = gen.sample_point(var.rng);
    p = p + c;
    std::list<Point> randPoints; //ds for storing rand points
    //use a large walk length e.g. 1000
    rand_point_generator(P, p, 1, 50*n, randPoints, var);
    // 3. Sample points from P
    int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
    randPoints.clear();
    rand_point_generator(P, p, num_of_samples, walk_len, randPoints, var);
    NT current_dist, max_dist;
    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist){
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    NT R=max_dist/radius;
    int mm=randPoints.size();
    boost::numeric::ublas::matrix<double> Ap(n,mm);
    for(int j=0; j<mm; j++){
        Point temp=randPoints.front();
        randPoints.pop_front();
        for (int i=0; i<n; i++){
            Ap(i,j)=double(temp[i]);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n);
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    double elleps=Minim::KhachiyanAlgo(Ap,0.01,w,Q,c2);

    Eigen::MatrixXd E(n,n);
    Eigen::VectorXd e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(int i=0; i<n; i++){
        e(i)=c2(i);
        for (int j=0; j<n; j++){
            E(i,j)=Q(i,j);
        }
    }

    int m=P.num_of_hyperplanes();
    Eigen::MatrixXd A(m,n);
    Eigen::VectorXd b(m);
    for(int i=0; i<m; ++i){
        b(i) = P.get_coeff(i,0);
        for(int j=1; j<n+1; ++j){
            A(i,j-1) = P.get_coeff(i,j);
        }
    }

    //Find the smallest and the largest axes of the elliposoid
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(E);
    NT rel = std::real(eigensolver.eigenvalues()[0]);
    NT Rel = std::real(eigensolver.eigenvalues()[0]);
    for(int i=1; i<n; i++){
        if(std::real(eigensolver.eigenvalues()[i])<rel) rel=std::real(eigensolver.eigenvalues()[i]);
        if(std::real(eigensolver.eigenvalues()[i])>Rel) Rel=std::real(eigensolver.eigenvalues()[i]);

    }
    //std::cout<<rel<<" "<<Rel<<std::endl;

    Eigen::LLT<Eigen::MatrixXd> lltOfA(E); // compute the Cholesky decomposition of E
    Eigen::MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    b = b - A*e;

    Eigen::MatrixXd L_1 = L.inverse();
    A = A*(L_1.transpose());

    // Write changes (actually perform rounding) to the polytope!
    for(int i=0; i<m; ++i){
        P.put_coeff(i,0,b(i));
        for(int j=1; j<n+1; ++j){
            P.put_coeff(i,j,A(i,j-1));
        }
    }

    //return L_1.determinant();
    return std::pair<double,double> (L_1.determinant(),rel/Rel);
}


// -------- ROTATION ---------- //
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
}

// -------- ROTATION ---------- //
template <class T>
double rotating(T &P){

  bool print = true; 
  //if(print) std::cout<<"\nRotate..."<<std::endl;
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
  //std::cout<<A<<"\n"<<b<<std::endl;
  
  MT M = MT::Random(d,d);
  //std::cout << "Here is the matrix m:" << std::endl << M << std::endl;
  Eigen::JacobiSVD<MT> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  //std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
  //std::cout << "Its left singular vectors are the columns of the U matrix:" << std::endl << svd.matrixU() << std::endl;
  //std::cout << "Det of U matrix:" << std::endl << svd.matrixU().determinant() << std::endl;
  //std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
  
  A = A*svd.matrixU();
  
  //std::cout<<A<<"\n"<<b<<std::endl;
  
  // Write changes (actually perform rotation) to the polytope!
	for(int i=0; i<m; ++i){
		P.put_coeff(i,0,b(i));
		for(int j=1; j<d+1; ++j){
		  P.put_coeff(i,j,A(i,j-1));
		}
	}

	//return std::abs(svd.matrixU().determinant());
	return std::abs(svd.matrixU().inverse().determinant());
}

#endif
