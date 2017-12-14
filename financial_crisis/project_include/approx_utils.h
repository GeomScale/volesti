
typedef boost::mt19937 RNGType;
//typedef CGAL::Gmpq                  EXACT_NT;
//typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel2;
//typedef Kernel2::RT					RT;

//typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
//typedef CGAL::Gmpq 			EXACT_NT;
//typedef double NT;
typedef CGAL::Cartesian_d<double> 	      Kernel; 
typedef Kernel::Point_d Point_d;
typedef CGAL::Timer				  timer;
typedef Kernel::Hyperplane_d 		Plane_d;

Point_d get_NewP2(Point_d p, double lamda, int k);

using namespace Eigen;

double metro_sq(Point_d A){
	
	int i,dim=A.dimension();
	double res=0.0;
	
	for (i=0; i<dim; i++){
		res+=std::pow(A[i],2);
	}
	
	return res;
	
}

void ray_facets(Point_d p0, int k, std::vector<double> &lamdas){
	
	double sum=0.0;
	int i,dim=p0.dimension();
	lamdas.push_back(-p0[k]);
	
	for (i=0; i<dim; i++){
		sum+=p0[i];
	}
	lamdas.push_back(1.0-sum);
	
	return;
	
}


void ray_NewFacets(Point_d p0, int k, std::vector<double> facet, double z1, double z2, std::vector<double> &lamdas){
	
	int i,dim=facet.size();
	double sum=0.0, lamda1, lamda2;
	
	for (i=0; i<dim; i++){
		sum+=facet[i]*p0[i];
	}
	
	lamda1=(z1-sum)/facet[k];
	lamda2=(z2-sum)/facet[k];
	
	lamdas.push_back(lamda1);
	lamdas.push_back(lamda2);
	
	return;	
	
}

std::pair<double,double> sect_ray_ball2(Point_d p0, Point_d center, double radius, int k){
	double a=1.0,b,c,D;
	int i,j,dim=p0.dimension();
	std::vector<double> temp_p;
	std::pair<double,double> result;
	Point_d p;
	
	for (i=0; i<dim; i++){
		temp_p.push_back(p0[i]-center[i]);
	}
	p=Point_d(dim, temp_p.begin(), temp_p.end());
	
	c=metro_sq(p)-std::pow(radius,2);
	b=2*p[k];
	
	Matrix2f A;// b_eig=b/a; c_eig=c/a;
	A(0,0)= 1.0; A(0,1)=1.0; A(1,0)=-(1.0+b+c); A(1,1)=-(1.0+b);
	EigenSolver<Matrix2f> es(A);
	double eig1 = es.eigenvalues()[0].real();
	double eig2 = es.eigenvalues()[1].real();
	//std::cout<<"[ball] eig1 is; "<<eig1<<" eig2 is: "<<eig2<<std::endl;
	
	D=std::pow(b,2)-4*a*c;
	//std::cout<<D<<" "<< (-b+std::sqrt(D))/(2*a)<<" "<<(-b-std::sqrt(D))/(2*a)<<std::endl;
	//lamdas.push_back((-b+std::sqrt(D))/(2*a));
	//lamdas.push_back((-b-std::sqrt(D))/(2*a));
	//lamdas.push_back(eig1);
	//lamdas.push_back(eig2);
	//result.first=eig2;
	//result.second=eig1;
	result.first=eig2;
	result.second=eig1;
	
	return result;
}


class ellipsoids{
	
	public:
		int dim;
		double c0;
		VectorXd center;
		MatrixXd C;
		
		ellipsoids(MatrixXd A, VectorXd b, double c, int d){
			dim=d;
			C=A;
			center=b;
			c0=c;
		}
		
		bool IsIn(Point_d p0){
			int i,dim=p0.dimension();
			VectorXd p(dim);
			for (i=0; i<dim; i++){
				p(i)=p0[i];
			}
			//std::cout<<(p-center).transpose()*C*(p-center)<<" "<<c0<<std::endl;
			if ((p-center).transpose()*C*(p-center)<c0){
				return true;
			}else{
				return false;
			}
		}
		
		void SectsRay(Point_d p, int k, std::vector<double> &lamdas){
			//SelfAdjointEigenSolver<Matrix2f> es;
			double a=0.0,b=0.0,c=0.0,D, b_eig, c_eig;
			int i,j,dim=p.dimension();
			
			for (i=0; i<dim; i++){
				if (i==k){
					b+=2*C(i,i)*(p[i]-center(i));
					a+=C(i,i);
				}
				c+=C(i,i)*std::pow(p[i]-center(i),2);
				for (j=i+1; j<dim; j++){
					if (i==k){
						b+=2*C(i,j)*(p[j]-center(j));
					}
					if (j==k){
						b+=2*C(i,j)*(p[i]-center(i));
					}
					c+=C(i,j)*(p[i]-center(i))*(p[j]-center(j))*2;
				}
			}
			c-=c0;
			//std::cout<<"k is: "<<k<<" a is: "<<a<<" b: "<<b<<" c: "<<c<<std::endl;
			//Matrix2f A; b_eig=b/a; c_eig=c/a;
			//A(0,0)= 1.0; A(0,1)=1.0; A(1,0)=-(1.0+b_eig+c_eig); A(1,1)=-(1.0+b_eig);
			//EigenSolver<Matrix2f> es(A);
		//	double eig1 = es.eigenvalues()[0].real();
		//	double eig2 = es.eigenvalues()[1].real();
			//std::cout<<"eig1 is; "<<eig1<<" eig2 is: "<<eig2<<std::endl;
			
			D=std::pow(b,2)-4*a*c;
			lamdas.push_back((-b+std::sqrt(D))/(2*a));
			lamdas.push_back((-b-std::sqrt(D))/(2*a));
			//lamdas.push_back(eig1);
			//lamdas.push_back(eig2);
			//std::cout<< (-b+std::sqrt(D))/(2*a)<<" "<<(-b-std::sqrt(D))/(2*a)<<std::endl;
			
			
			return;
		}
		
		
		void SectsRay2(Point_d p, int k, std::vector<double> &lamdas){
			//SelfAdjointEigenSolver<Matrix2f> es;
			double a=0.0,b=0.0,c=0.0,D, b_eig, c_eig, x1,x2;
			int i,j,dim=p.dimension();
			
			for (i=0; i<dim; i++){
				if (i==k){
					b+=2*C(i,i)*(p[i]-center(i));
					a+=C(i,i);
				}
				c+=C(i,i)*std::pow(p[i]-center(i),2);
				for (j=i+1; j<dim; j++){
					if (i==k){
						b+=2*C(i,j)*(p[j]-center(j));
					}
					if (j==k){
						b+=2*C(i,j)*(p[i]-center(i));
					}
					c+=C(i,j)*(p[i]-center(i))*(p[j]-center(j))*2;
				}
			}
			c-=c0;
			//std::cout<<"k is: "<<k<<" a is: "<<a<<" b: "<<b<<" c: "<<c<<std::endl;
			//Matrix2f A; b_eig=b/a; c_eig=c/a;
			//A(0,0)= 1.0; A(0,1)=1.0; A(1,0)=-(1.0+b_eig+c_eig); A(1,1)=-(1.0+b_eig);
			//EigenSolver<Matrix2f> es(A);
		//	double eig1 = es.eigenvalues()[0].real();
		//	double eig2 = es.eigenvalues()[1].real();
			//std::cout<<"eig1 is; "<<eig1<<" eig2 is: "<<eig2<<std::endl;
			
			D=std::pow(b,2)-4*a*c;
			//std::cout<<D<<" "<< (-b+std::sqrt(D))/(2*a)<<" "<<(-b-std::sqrt(D))/(2*a)<<" "<<IsIn(p)<<std::endl;
			if (D>=0.0){
				
				x1=(-b+std::sqrt(D))/(2*a);
				x2=(-b-std::sqrt(D))/(2*a);
				lamdas.push_back(x1);
				lamdas.push_back(x2);
			}
			//std::cout<<"x1 is; "<<x1<<" x2 is: "<<x2<<std::endl;
			/*if (x1<0.0 && x2<0.0){
				if (x1>x2){
					lamdas.push_back(x1);
				}else{
					lamdas.push_back(x2);
				}
			}else{
				if (x1<x2){
					lamdas.push_back(x1);
				}else{
					lamdas.push_back(x2);
				}
			}
			}*/
			//std::cout<<"x1 is; "<<x1<<" x2 is: "<<x2<<std::endl;
			//if (D>=0.0){
			//	lamdas.push_back((-b+std::sqrt(D))/(2*a));
			//	lamdas.push_back((-b-std::sqrt(D))/(2*a));
			//}
			//lamdas.push_back(eig1);
			//lamdas.push_back(eig2);
			//std::cout<< (-b+std::sqrt(D))/(2*a)<<" "<<(-b-std::sqrt(D))/(2*a)<<std::endl;
			
			
			return;
		}
		
		
		std::vector<double> getConst(){
			int i;
			double max=0.0,step,min=0.0;
			min=0.6843;
			max=0.6905;
			std::vector<double> consts(99);
			//for (i=0; i<dim; i++){
			//	if (C(i,i)>max){
				//	max=C(i,i);
				//}
				//if (C(i,i)<min){
					//min=C(i,i);
				//}
				
			//}
			//std::cout<<"min max"<<min<<" "<<max::std::endl;
			step=(max-min)/100.0;
			for (i=1; i<=99; i++){
				consts[i-1]=min+i*step;
				//std::cout<<consts[i-1]<<" ";
			}
			return consts;
			
		}
		
		double MatMult(Point_d p0){
			int i,dim=p0.dimension();
			VectorXd p(dim);
			for (i=0; i<dim; i++){
				p(i)=p0[i];
			}
			//std::cout<<(p-center).transpose()*C*(p-center)<<" "<<c0<<std::endl;
			return (p-center).transpose()*C*(p-center);
			
		}
		
		int getDim(){
			return dim;
		}
		
		MatrixXd getMat(){
			return C;
		}
		
		VectorXd getCenter(){
			return center;
		}
	
};


Point_d get_NewP2(Point_d p, double lamda, int k){
	Point_d pnew;
	std::vector<double> temp_p;
	int dim=p.dimension(),i;
	
	for (i=0; i<dim; i++){
		if (i==k){
			temp_p.push_back(p[i]+lamda);
		}else{
			temp_p.push_back(p[i]);
		}
	}
	pnew=Point_d(dim, temp_p.begin(), temp_p.end());
	
	return pnew;
	
}


Point_d Subtrack_points2(Point_d p1, Point_d p2){
	int i,dim=p1.dimension();
	std::vector<double> temp_p;
	Point_d pnew;
	
	for(i=0; i<dim; i++){
		temp_p.push_back(p1[i]-p2[i]);
	}
	pnew=Point_d(dim, temp_p.begin(), temp_p.end());
	return pnew;
	
}


Point_d Add_points2(Point_d p1, Point_d p2){
	int i,dim=p1.dimension();
	std::vector<double> temp_p;
	Point_d pnew;
	
	for(i=0; i<dim; i++){
		temp_p.push_back(p1[i]+p2[i]);
	}
	pnew=Point_d(dim, temp_p.begin(), temp_p.end());
	return pnew;
	
}


bool isin_ball2(Point_d center, double radius, Point_d p){
	
	Point_d pnew;
	pnew=Subtrack_points2(p,center);
	if ( std::sqrt(metro_sq(pnew))<radius ){
		return true;
	}else{
		return false;
	}
}


std::pair<Point_d,double> get_center_radius_inscribed_simplex(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end){
	
	Point_d p0=*it_beg,p1,c;
	int dim=p0.dimension(),i,j;
	std::vector<double> temp_p;
	double radius=0.0,gi,sum=0.0;
	MatrixXd B(dim,dim);
	MatrixXd Bg(dim,dim);
	VectorXd e(dim);
	VectorXd row(dim);
	VectorXd g(dim);
	std::pair<Point_d,double> result;
	//MatrixXd res2(1,4);
	
	
	for (j=1; j<dim+1; j++){
		Point_d pk=*(it_beg+j);
		e(j-1)=1.0;
		for (i=0; i<dim; i++){
				B(i,j-1)=pk[i]-p0[i];
		}
	}
	Bg=B;
	B=B.inverse();
	//std::cout<< B <<std::endl;
	//std::cout<<"\n";
	for (i=0; i<dim; i++){
		for (j=0; j<dim; j++){
			row(j)=B(i,j);
		}
	//	std::cout<< row <<std::endl;
	//	std::cout<<"\n";
		gi=row.norm();
		radius+=gi;
		g(i)=gi;
		if (i<dim-1){
			sum+=gi;
		}
	//	std::cout<<radius<<std::endl;
	//	std::cout<<"\n";
	}
//	std::cout<<"hello"<<std::endl;
	e=e*B;
//	std::cout<< e <<std::endl;
	//std::cout<<"\n";
	radius+=e.norm();
	//std::cout<<"hello2.5"<<std::endl;
	radius=1.0/radius;
//	std::cout<<"hello2"<<std::endl;
	g=Bg*g;
	g=radius*g;
	//std::cout<<"hello3"<<std::endl;
	for (i=0; i<dim; i++){
		temp_p.push_back(p0[i]+g(i));
	}
	c=Point_d(dim,temp_p.begin(), temp_p.end());
	//std::cout<< "small radius is: "<< radius <<std::endl;
	result.first=c;
	result.second=radius;
	//std::cout<<"hello2"<<std::endl;
	return result;
	
}




std::pair<Point_d,double> rand_inscribed_ball(ellipsoids G){//std::vector<Point_d> &points, ellipsoids G){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,pointer,counter=0,dim=G.getDim();  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<Point_d> points;
	Point_d temp,center;
	bool inside;
	double radius;
	std::pair<Point_d,double> result;
	//std::vector<bool> filter;
	
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	rng.seed(static_cast<unsigned int>(std::time(0)));

	

	if (dim<=60){
		
		x_vec.assign(dim+1,0);

		while (counter<dim+1){

			// Set all the point's coordinates equal to zero and clear the integers' list
			y.assign(dim,0); pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){

				x_rand = uidist(rng);
				// Check if this integer is selected first time
				if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand) == x_vec.begin()+pointer ){
					pointer++;
					x_vec[pointer]=x_rand;
				}

			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );

			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			temp=Point_d(dim,y.begin(),y.end());
			inside=G.IsIn(temp);
			if (inside){
				counter++;
				points.push_back(temp);
			}
		
		}
	}else if(dim<=80){	
		
		x_vec.assign(dim+1,0);

		std::vector<bool> filter;
		bool t1=true,t2=true;
		pr=3*dim+1;
		if(pr%2==0) pr+=1;

		while(t1){
			t2=true;
			divisors=(int)floor(sqrt((double)pr))+1;
			for (i=3; i<divisors+1; i+=2){
				if (pr%i==0){
					t2=false;
					break;
				}
			}
			if(t2) break;
			pr+=2;
		}
		
		// Generate the number of points requested
		while(counter<dim+1){

			// Set all the point's coordinates equal to zero and the filter equal to true
			y.assign(dim,0); filter.assign(pr,true);
			pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){
				x_rand = uidist(rng);

				// Check if this integer is the first that is mapped to the specific filter's position
				if ( filter[x_rand%pr] ){
					filter[x_rand%pr]=false; pointer++;
					x_vec[pointer]=x_rand;
				}
			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );
	
			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			temp=Point_d(dim,y.begin(),y.end());
			inside=G.IsIn(temp);
			if (inside){
				counter++;
				points.push_back(temp);
			}
		
		}
	}else{

		std::vector<double> x_vec2;
		double Ti,sum;	

		x_vec2.assign(dim+1,0.0);
	
		// Generate the number of points requested
		while(counter<dim+1){

			// Set all the point's coordinates equal to zero and clear the integers' list
			pointer=0; sum=0.0;
			
			while ( pointer<dim+1 ){
	
				x_rand = uidist(rng);			
				Ti=-log(((double)x_rand)/((double)M));
				sum+=Ti;
				x_vec2[pointer]= Ti; pointer++;
	
			}
			
			for (j=0; j<dim; j++){
				x_vec2[j]/=sum;
			}
	
			points.push_back(Point_d(dim,x_vec2.begin(),x_vec2.end()-1));
			temp=Point_d(dim,x_vec2.begin(),x_vec2.end()-1);
			inside=G.IsIn(temp);
			if (inside){
				counter++;
				points.push_back(temp);
			}
			
		}

	}
	std::cout<<"size of points: "<<points.size()<<std::endl;
	result=get_center_radius_inscribed_simplex(points.begin(), points.end());
	std::cout<<"center is: "<<std::endl;
	//std::cout<<"radius is: "<<*r<<std::endl;
	//*radius=*r;
	//c=center;

	return result;

}




std::pair<Point_d,double> rand_inscribed_ball2(ellipsoids G, std::vector<double> plane, double z1, double z2){//std::vector<Point_d> &points, ellipsoids G){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,pointer,counter=0,dim=G.getDim();  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<Point_d> points;
	Point_d temp,center;
	bool inside;
	double radius;
	std::pair<Point_d,double> result;
	//std::vector<bool> filter;
	
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	rng.seed(static_cast<unsigned int>(std::time(0)));

	

	if (dim<=60){
		double sumin;
		x_vec.assign(dim+1,0);

		while (counter<dim+1){

			// Set all the point's coordinates equal to zero and clear the integers' list
			y.assign(dim,0); pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){

				x_rand = uidist(rng);
				// Check if this integer is selected first time
				if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand) == x_vec.begin()+pointer ){
					pointer++;
					x_vec[pointer]=x_rand;
				}

			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );

			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			temp=Point_d(dim,y.begin(),y.end());
			inside=G.IsIn(temp);
			if (inside){
				sumin=0.0;
				for (j=0; j<dim; j++){
					sumin+=plane[j]*temp[j];
				}
				if (sumin<z2 && sumin>z1){
					counter++;
					points.push_back(temp);
				}
			}
		
		}
	}else if(dim<=80){	
		
		x_vec.assign(dim+1,0);

		std::vector<bool> filter;
		bool t1=true,t2=true;
		pr=3*dim+1;
		if(pr%2==0) pr+=1;

		while(t1){
			t2=true;
			divisors=(int)floor(sqrt((double)pr))+1;
			for (i=3; i<divisors+1; i+=2){
				if (pr%i==0){
					t2=false;
					break;
				}
			}
			if(t2) break;
			pr+=2;
		}
		
		// Generate the number of points requested
		while(counter<dim+1){

			// Set all the point's coordinates equal to zero and the filter equal to true
			y.assign(dim,0); filter.assign(pr,true);
			pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){
				x_rand = uidist(rng);

				// Check if this integer is the first that is mapped to the specific filter's position
				if ( filter[x_rand%pr] ){
					filter[x_rand%pr]=false; pointer++;
					x_vec[pointer]=x_rand;
				}
			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );
	
			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			temp=Point_d(dim,y.begin(),y.end());
			inside=G.IsIn(temp);
			if (inside){
				counter++;
				points.push_back(temp);
			}
		
		}
	}else{

		std::vector<double> x_vec2;
		double Ti,sum;	

		x_vec2.assign(dim+1,0.0);
	
		// Generate the number of points requested
		while(counter<dim+1){

			// Set all the point's coordinates equal to zero and clear the integers' list
			pointer=0; sum=0.0;
			
			while ( pointer<dim+1 ){
	
				x_rand = uidist(rng);			
				Ti=-log(((double)x_rand)/((double)M));
				sum+=Ti;
				x_vec2[pointer]= Ti; pointer++;
	
			}
			
			for (j=0; j<dim; j++){
				x_vec2[j]/=sum;
			}
	
			points.push_back(Point_d(dim,x_vec2.begin(),x_vec2.end()-1));
			temp=Point_d(dim,x_vec2.begin(),x_vec2.end()-1);
			inside=G.IsIn(temp);
			if (inside){
				counter++;
				points.push_back(temp);
			}
			
		}

	}
	std::cout<<"size of points: "<<points.size()<<std::endl;
	result=get_center_radius_inscribed_simplex(points.begin(), points.end());
	//std::cout<<"center is: "<<center<<std::endl;
	//std::cout<<"radius is: "<<*r<<std::endl;
	//*radius=*r;
	//c=center;

	return result;

}









std::pair<Point_d,double> rand_inscribed_ball_nonConv(ellipsoids G1, ellipsoids G2, std::vector<double> plane, double z1, double z2){//std::vector<Point_d> &points, ellipsoids G){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,pointer,counter=0,dim=G1.getDim();  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<Point_d> points;
	Point_d temp,center;
	bool inside1,inside2;
	double radius, sumin;
	std::pair<Point_d,double> result;
	//std::vector<bool> filter;
	
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	rng.seed(static_cast<unsigned int>(std::time(0)));

	

	if (dim<=60){
		
		x_vec.assign(dim+1,0);

		while (counter<dim+1){

			// Set all the point's coordinates equal to zero and clear the integers' list
			y.assign(dim,0); pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){

				x_rand = uidist(rng);
				// Check if this integer is selected first time
				if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand) == x_vec.begin()+pointer ){
					pointer++;
					x_vec[pointer]=x_rand;
				}

			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );

			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			temp=Point_d(dim,y.begin(),y.end());
			inside1=G1.IsIn(temp);
			inside2=G2.IsIn(temp);
			
			if (inside1 && !inside2){
				sumin=0.0;
				for (j=0; j<dim; j++){
					sumin+=plane[j]*temp[j];
				}
				if (sumin<z2 && sumin>z1){
					counter++;
					points.push_back(temp);
					center=temp;
					radius=0.000000000001;
					result.first=center;
					result.second=radius;
					return result;
				}
			}
		
		}
	}else if(dim<=80){	
		
		x_vec.assign(dim+1,0);

		std::vector<bool> filter;
		bool t1=true,t2=true;
		pr=3*dim+1;
		if(pr%2==0) pr+=1;

		while(t1){
			t2=true;
			divisors=(int)floor(sqrt((double)pr))+1;
			for (i=3; i<divisors+1; i+=2){
				if (pr%i==0){
					t2=false;
					break;
				}
			}
			if(t2) break;
			pr+=2;
		}
		
		// Generate the number of points requested
		while(counter<dim+1){

			// Set all the point's coordinates equal to zero and the filter equal to true
			y.assign(dim,0); filter.assign(pr,true);
			pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){
				x_rand = uidist(rng);

				// Check if this integer is the first that is mapped to the specific filter's position
				if ( filter[x_rand%pr] ){
					filter[x_rand%pr]=false; pointer++;
					x_vec[pointer]=x_rand;
				}
			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );
	
			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			temp=Point_d(dim,y.begin(),y.end());
			inside1=G1.IsIn(temp);
			inside2=G2.IsIn(temp);
			if (inside1 && !inside2){
				sumin=0.0;
				for (j=0; j<dim; j++){
					sumin+=plane[j]*temp[j];
				}
				if (sumin<z2 && sumin>z1){
					counter++;
					points.push_back(temp);
				}
			}
		
		}
	}else{

		std::vector<double> x_vec2;
		double Ti,sum;	

		x_vec2.assign(dim+1,0.0);
	
		// Generate the number of points requested
		while(counter<dim+1){

			// Set all the point's coordinates equal to zero and clear the integers' list
			pointer=0; sum=0.0;
			
			while ( pointer<dim+1 ){
	
				x_rand = uidist(rng);			
				Ti=-log(((double)x_rand)/((double)M));
				sum+=Ti;
				x_vec2[pointer]= Ti; pointer++;
	
			}
			
			for (j=0; j<dim; j++){
				x_vec2[j]/=sum;
			}
	
			points.push_back(Point_d(dim,x_vec2.begin(),x_vec2.end()-1));
			temp=Point_d(dim,x_vec2.begin(),x_vec2.end()-1);
			inside1=G1.IsIn(temp);
			inside2=G2.IsIn(temp);
			if (inside1 && !inside2){
				sumin=0.0;
				for (j=0; j<dim; j++){
					sumin+=plane[j]*temp[j];
				}
				if (sumin<z2 && sumin>z1){
					counter++;
					points.push_back(temp);
				}
			}
			
		}

	}
	std::cout<<"size of points: "<<points.size()<<std::endl;
	result=get_center_radius_inscribed_simplex(points.begin(), points.end());
	//std::cout<<"center is: "<<center<<std::endl;
	//std::cout<<"radius is: "<<*r<<std::endl;
	//*radius=*r;
	//c=center;

	return result;

}

