//Sample from the unit simplex

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

using namespace Eigen;

void Sam_Unit(int dim, int num, std::vector<Point_d> &points){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	//std::vector<bool> filter;
	
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	rng.seed(static_cast<unsigned int>(std::time(0)));

	

	if (dim<=60){
		
		x_vec.assign(dim+1,0);

		for (i=0; i<num; i++){

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
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
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
		for (i=0; i<num; i++){

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
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}else{

		std::vector<double> x_vec2;
		double Ti,sum;	

		x_vec2.assign(dim+1,0.0);
	
		// Generate the number of points requested
		for (i=0; i<num; i++){

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
			
		}

	}

	return;

}


//sample from Facet simplex  //Needs adjustments
void Sam_Unit_NoProjection(int dim, int num, std::vector<Point_d> &points){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	//std::vector<bool> filter;
	dim--;
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	rng.seed(static_cast<unsigned int>(std::time(0)));

	

	/*if (dim<=60){
		
		x_vec.assign(dim+1,0);

		for (i=0; i<num; i++){

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
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
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
		for (i=0; i<num; i++){

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
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}else{*/

		std::vector<double> x_vec2;
		double Ti,sum;	

		x_vec2.assign(dim+1,0.0);
	
		// Generate the number of points requested
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and clear the integers' list
			pointer=0; sum=0.0;
			
			while ( pointer<dim+1 ){
	
				x_rand = uidist(rng);			
				Ti=-log(((double)x_rand)/((double)M));
				sum+=Ti;
				x_vec2[pointer]= Ti; pointer++;
	
			}
			
			for (j=0; j<dim+1; j++){
				x_vec2[j]/=sum;
			}
	
			points.push_back(Point_d(dim+1,x_vec2.begin(),x_vec2.end()));
			
		}

	//}

	return;

}



//Owen mapping for sample from an arbitrary simplex given in V-represantation
void Sam_arbest(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end, int num, std::vector<Point_d> &points){
	
	int n=std::distance(it_beg,it_end),dim,j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	//std::vector<double> x_vec2;
	
	double Xj;
	

	Point_d p0=*it_beg;
	if (n>=2){   // Check if points could define a simplex
		if (p0.dimension()!=n-1){   // Check if the polytope is a Simplex
			std::cout<<"This Polytope is not a Simplex. The problem is open."<<std::endl;
			return;
		}else{
			dim=p0.dimension();
		}
	}else{
		std::cout<<"We need more points"<<std::endl;
		return;
	}

	
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	rng.seed(static_cast<unsigned int>(std::time(0)));

	

	if (dim<=60){
		
		x_vec.assign(dim+2,0);  x_vec[dim+1]=M;

		for (i=0; i<num; i++){

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
			for(j=0; j<dim+1; j++){
				Point_d pk=*(it_beg+j);
				Xj=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for (k=0; k<dim; k++){
					y[k]+=Xj*pk[k];
				}
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}else if(dim<=80){	

		//printf("hello_filter");

		std::vector<bool> filter;
		
		x_vec.assign(dim+2,0);  x_vec[dim+1]=M;
		
		//double Xj;
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
		for (i=0; i<num; i++){

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
			for(j=0; j<dim+1; j++){
				Point_d pk=*(it_beg+j);
				Xj=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for (k=0; k<dim; k++){
					y[k]+=Xj*pk[k];
				}
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}else{
		//printf("hello0.5");
		std::vector<double> x_vec2;
		double Ti,sum;	
		//printf("hello0.5");
		x_vec2.assign(dim+1,0.0);
		//printf("hello0.6");
	
		// Generate the number of points requested
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and clear the integers' list
			pointer=0; sum=0.0; y.assign(dim,0);
			//printf("hello1");
			while ( pointer<dim+1 ){
	
				x_rand = uidist(rng);			
				Ti=-log(((double)x_rand)/((double)M));
				sum+=Ti;
				x_vec2[pointer]= Ti; pointer++;
	
			}

			//printf("hello2");
			
			for (j=0; j<dim+1; j++){
				x_vec2[j]/=sum;
				Point_d pk=*(it_beg+j);
				Xj=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for (k=0; k<dim; k++){
					y[k]+=x_vec2[j]*pk[k];
				}
			}
			//printf("hello3");
	
			points.push_back(Point_d(dim,y.begin(),y.end()));
			
		}

	}

	return;

}


std::pair<double,double> hit_and_run_newP2(ellipsoids G,Point_d p0, Point_d center, double radius, int k, bool isball, bool &onBall){
	
	int i,dim=p0.dimension(),n;
	double lamda1, lamda2, lamdaball1, lamdaball2;
	std::pair<double,double> result, lamdasball;
	std::vector<double> temp_p, lamdas;
	
	/*if (!first){
		for (i=0; i<dim; i++){
			if (i!=kold){
				temp_p.push_back(p0[i]);
			}else{
				temp_p.push_back(p0[i]+lamda);
			}
		}
		p0=Point_d(dim, temp_p.begin(), temp_p.end());
	}*/
	
	ray_facets(p0, k, lamdas);
	G.SectsRay(p0, k, lamdas);
	//ray_NewFacets(p0, k, facet1, z1, z2, lamdas);
	//if (isball){
	//	sect_ray_ball(p0, center, radius, k, lamdas);
	//}
	
	n=lamdas.size();
	lamda1=-1234976.34;
	lamda2=12384.34;
	
	for (i=0; i<n; i++){
		if (lamdas[i]<0){
			if (lamdas[i]>lamda1){
				lamda1=lamdas[i];
			}
		}else{
			if (lamdas[i]<lamda2){
				lamda2=lamdas[i];
			}
		}
	}
	
	if (isball){
		lamdasball=sect_ray_ball2(p0, center, radius, k);
		lamdaball1=lamdasball.first;
		lamdaball2=lamdasball.second;
		
		if(lamdaball1>lamda1){
			lamda1=lamdaball1;
			onBall=(onBall && true);
		}else{
			onBall=false;
			//std::cout<<"make it false\n";
		}
		if(lamdaball2<lamda2){
			lamda2=lamdaball2;
			onBall=(onBall && true);
		}else{
			onBall=false;
			//std::cout<<"make it false\n";
		}
	}
	
	result.first=lamda1;
	result.second=lamda2;
	
	return result;
}



//Hit and run for unit simplex and  one ellispoid intersection
std::pair<double,double> hit_and_run_newP3(ellipsoids G,Point_d p0, Point_d center, double radius, std::vector<double> facet, double z1, double z2, int k, bool isball, bool &onBall){
	
	int i,dim=p0.dimension(),n;
	double lamda1, lamda2, lamdaball1, lamdaball2;
	std::pair<double,double> result, lamdasball;
	std::vector<double> temp_p, lamdas;
	
	/*if (!first){
		for (i=0; i<dim; i++){
			if (i!=kold){
				temp_p.push_back(p0[i]);
			}else{
				temp_p.push_back(p0[i]+lamda);
			}
		}
		p0=Point_d(dim, temp_p.begin(), temp_p.end());
	}*/
	
	ray_facets(p0, k, lamdas);
	G.SectsRay(p0, k, lamdas);
	ray_NewFacets(p0, k, facet, z1, z2, lamdas);
	//if (isball){
	//	sect_ray_ball(p0, center, radius, k, lamdas);
	//}
	
	n=lamdas.size();
	lamda1=-1234976.34;
	lamda2=12384.34;
	
	for (i=0; i<n; i++){
		if (lamdas[i]<0){
			if (lamdas[i]>lamda1){
				lamda1=lamdas[i];
			}
		}else{
			if (lamdas[i]<lamda2){
				lamda2=lamdas[i];
			}
		}
	}
	
	if (isball){
		lamdasball=sect_ray_ball2(p0, center, radius, k);
		lamdaball1=lamdasball.first;
		lamdaball2=lamdasball.second;
		
		if(lamdaball1>lamda1){
			lamda1=lamdaball1;
			onBall=(onBall && true);
		}else{
			onBall=false;
			//std::cout<<"make it false\n";
		}
		if(lamdaball2<lamda2){
			lamda2=lamdaball2;
			onBall=(onBall && true);
		}else{
			onBall=false;
			//std::cout<<"make it false\n";
		}
	}
	
	result.first=lamda1;
	result.second=lamda2;
	
	return result;
}


//hit and run for nonconvex body defined by unit siimplex, two ellipoids and two parallel hyperplanes
std::pair<double,double> hit_and_run_newP_nonConv(ellipsoids G1, ellipsoids G2, Point_d p0, Point_d center, double radius, std::vector<double> facet, double z1, double z2, int k, bool isball, bool &onBall){
	
	int i,dim=p0.dimension(),n;
	double lamda1, lamda2, lamdaball1, lamdaball2;
	std::pair<double,double> result, lamdasball;
	std::vector<double> temp_p, lamdas;
	
	/*if (!first){
		for (i=0; i<dim; i++){
			if (i!=kold){
				temp_p.push_back(p0[i]);
			}else{
				temp_p.push_back(p0[i]+lamda);
			}
		}
		p0=Point_d(dim, temp_p.begin(), temp_p.end());
	}*/
	
	ray_facets(p0, k, lamdas);
	G1.SectsRay(p0, k, lamdas);
	ray_NewFacets(p0, k, facet, z1, z2, lamdas);
	G2.SectsRay2(p0, k, lamdas);
	//if (isball){
	//	sect_ray_ball(p0, center, radius, k, lamdas);
	//}
	
	n=lamdas.size();
	lamda1=-1234976.34;
	lamda2=12384.34;
	
	for (i=0; i<n; i++){
		if (lamdas[i]<0){
			if (lamdas[i]>lamda1){
				lamda1=lamdas[i];
			}
		}else{
			if (lamdas[i]<lamda2){
				lamda2=lamdas[i];
			}
		}
	}
	
	if (isball){
		lamdasball=sect_ray_ball2(p0, center, radius, k);
		lamdaball1=lamdasball.first;
		lamdaball2=lamdasball.second;
		
		if(lamdaball1>lamda1){
			lamda1=lamdaball1;
			onBall=(onBall && true);
		}else{
			onBall=false;
			//std::cout<<"make it false\n";
		}
		if(lamdaball2<lamda2){
			lamda2=lamdaball2;
			onBall=(onBall && true);
		}else{
			onBall=false;
			//std::cout<<"make it false\n";
		}
	}
	
	result.first=lamda1;
	result.second=lamda2;
	
	return result;
}

