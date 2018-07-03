// VolEsti

// Copyright (c) 2012-2017 Vissarion Fisikopoulos

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

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
    return std::pair<Point,NT> (c,R);
}


template <class T1>
NT rounding(T1 &P , Point c, NT radius, vars &var){
    //typedef typename T1::FT 	K;
    int n=var.n, walk_len=var.walk_steps;
    // 1. Compute the Chebychev ball (largest inscribed ball) with center and radius 
	//Point c(n);       //center
    //K radius;
    //P.chebyshev_center(c,radius);
    
    // 2. Generate the first random point in P
  // Perform random walk on random point in the Chebychev ball 
	//if (print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
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
    
    
    // 4. Compute the transformation matrix T
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(n,n);
    bool well_rounded=false;
    //int t=8*n*n*n,
    int t=var.m;
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
        rand_point_generator(P, p, 1, 50*n, randPoints, var);
        randPoints.clear();
        rand_point_generator(P2, p, t, walk_len, randPoints, var);
        Eigen::MatrixXd PM=getPointsMat(randPoints,n);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(PM, Eigen::ComputeThinU | Eigen::ComputeThinV);

        NT min=svd.singularValues()(0);
        for(int i=1; i<n; i++){
            if(svd.singularValues()(i)<min){
                min=svd.singularValues()(i);
            }
        }
        for(int i=0; i<n; i++){
            S(i,i)=svd.singularValues()(i)/min;
        }

        T=svd.matrixV()*S.inverse()*T;
        well_rounded=true;
        for (int i=0; i<n; i++){
            if (svd.singularValues()(i)>2.0){
                if (tries>(int)std::log(R)){
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
    
    return T.inverse().determinant();
}


// -------- ROTATION ---------- //
template <class T>
double rotating_old(T &P){
    
    return 1.0;
}
/*
	
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
template <class T>
double rotating(T &P){
    
    return 1.0;
}
/*
	
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
  
	return 0;
}*/


// ----- ROUNDING ------ //
template <class T>
double rounding(T &P, 
                vars &var,  // constans for volume
                vars &var2 // constants for optimization in case of MinkSums
                ){
    return 1.0;
}
/*
  bool print = var.verbose;
  int n = var.n;
  int rnum = var.m;
	int walk_len = var.walk_steps;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> urdist = var.urdist;
	boost::random::uniform_int_distribution<> uidist(0,n-1);

  // 1. Compute the Chebychev ball (largest inscribed ball) with center and radius 
	Point c;       //center
    double radius;
    P.chebyshev_center(c,radius);
    if (print) std::cout<<"Chebychev center= "<<c<<"\nradius="<<radius<<std::endl;
  
  // 2. Generate the first random point in P
  // Perform random walk on random point in the Chebychev ball 
	if (print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
	CGAL::Random_points_in_ball_d<Point> gen (n, radius);
	Point p = *gen;
	p = p + (c-CGAL::Origin());
	std::list<Point> randPoints; //ds for storing rand points
	//use a large walk length e.g. 1000
	rand_point_generator(P, p, 1, 1000, randPoints, var); 
	if (print) std::cout<<"First random point: "<<p<<std::endl;
		
	// 3. Sample points from P
	//randPoints.push_front(p);
	int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
	if(print) std::cout<<"\nCompute "<<num_of_samples<<" random points in P"<<std::endl;
	rand_point_generator(P, p, num_of_samples, walk_len, randPoints, var); 
	
    // 4. Compute approximation of min enclosing ellipsoid of randPoints
    if(print) std::cout<<"\nCompute approximate min ellipsoid..."<<std::endl;
    Traits traits;
    AME ame(0.01, randPoints.begin(), randPoints.end(), traits);
	//std::cout<<ame.defining_matrix(1,1)<<std::endl;
	typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VT;
    // Construct polytope matrices
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
	// 4a. Construct ellipsoid matrix
	int k=randPoints.size();
	double achieved_epsilon = ame.achieved_epsilon();
	MT E(d,d);
	for(int i=0; i<d; ++i){
		for(int j=0; j<d; ++j){
		  E(i,j) = ame.defining_matrix(i,j) / ((1+achieved_epsilon)*(1+d));
		}
	}
	//std::cout<<E<<std::endl;
	// 4b. The center of the ellipsoid
	VT c_e(d);
	for(AME::Center_coordinate_iterator cit=ame.center_cartesian_begin(); cit!=ame.center_cartesian_end(); ++cit)
		c_e(cit-ame.center_cartesian_begin()) = *cit;
	//std::cout<<"center\n"<<c_e<<std::endl;
	
	// The defining vector of the ellipsoid
	VT e(d);
	for(int i=0; i<d; ++i)
		e(i) = ame.defining_vector(i) / ((1+achieved_epsilon)*(1+d));
	//std::cout<<"e:\n"<<e<<std::endl;
	
	//std::cout<<"test e:\n"<<(-1*E.transpose()*c_e)-E*c_e<<std::endl;
	
	//std::cout<<"Is full dimensional: "<<ame.is_full_dimensional()<<std::endl;
	// Axes lengths
	/*
	AME::Axes_lengths_iterator axes = ame.axes_lengths_begin();
	for (int i = 0; i < d; ++i) {
		std::cout << "Semiaxis " << i << " has length " << *axes++  << "\n"
							<< "and Cartesian coordinates ";
		for (AME::Axes_direction_coordinate_iterator
					 d_it = ame.axis_direction_cartesian_begin(i);
				 d_it != ame.axis_direction_cartesian_end(i); ++d_it)
			std::cout << *d_it << ' ';
		std::cout << ".\n";
	}
	*/
    /*
	Eigen::LLT<MT> lltOfA(E);  compute the Cholesky decomposition of E
	MT L = lltOfA.matrixL();  retrieve factor L  in the decomposition
	//std::cout<<L<<std::endl;
	
	//std::cout<<L*L.transpose()<<std::endl;
	
	//b = b - A*c_e;
	
	//MT L_1 = L.inverse();
	A = A*(L_1.transpose());
	
	// Write changes (actually perform rounding) to the polytope!
	for(int i=0; i<m; ++i){
		P.put_coeff(i,0,b(i));
		for(int j=1; j<d+1; ++j){
		  P.put_coeff(i,j,A(i,j-1));
		}
	}
	
	//const double pi = boost::math::constants::pi<double>();
	//double ame_vol = (std::pow(pi,d/2.0) / (std::tgamma((d/2.0)+1) * std::sqrt(E.determinant())));
	//std::cout<<"Ellipsoid volume "<<ame_vol<<" ,"<<std::pow(2,d)<<std::endl;
	//std::cout<<"rounding value="<<L_1.determinant()<<" , "<<L_1.determinant()*std::pow(2,d)<<std::endl;	
	
	return L_1.determinant();
}*/
