typedef CGAL::Gmpq                  EXACT_NT;

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
//double VolEsti_ellips2(ellipsoids G, int walk_len, int rnum);


// VolEsti for computing 
double VolEsti_ellips2(ellipsoids G, int walk_len, int rnum){
	
	double radius,lamda1,lamda2,lamda,l_rand,radLarge,radSmall;
	int dim=G.getDim(),i,j,k;
	Point_d c;       //center
	std::pair<double,double> lamdas;
	std::pair<Point_d,double> inscribed_ball;
	bool ball_inside=false, ray_onball=true;
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// the random engine with this seed
	RNGType rng(seed);
	//RNGType rng;
	boost::random::uniform_int_distribution<> uidist1(0,dim-1);
	boost::random::uniform_real_distribution<> urdist2(0,1); 
	
	//get the inscribed ball
	inscribed_ball=rand_inscribed_ball(G);
	c=inscribed_ball.first;
	radius=inscribed_ball.second;
	std::cout<<"center is: "<<c<<std::endl;
	
	//radius=*rad;
	std::cout<<"radius is: "<<radius<<std::endl;
	
	CGAL::Random_points_in_ball_d<Point_d> gen (dim, radius);
    Point_d p = *gen;
    //p = p + (c-CGAL::Origin());
    p=Add_points2(p,c);
    std::vector<Point_d> randPoints;
    std::cout<<"first rand point to begin HnR in P: "<<p<<std::endl;
    std::cout<<G.IsIn(p)<<std::endl;
	
	k = uidist1(rng); //rand coord
	
	//for(i=0; i<
	lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball); //needs changes
	lamda1=lamdas.first;
	lamda2=lamdas.second;
	l_rand=urdist2(rng);
	lamda=lamda1+(lamda2-lamda1)*l_rand;
	p=get_NewP2(p, lamda, k);
	std::cout<<"second rand point to begin HnR in P: "<<p<<std::endl;
	for (i=0; i<1000; i++){
		k = uidist1(rng);
		std::cout<<k<<std::endl;
		lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
		lamda1=lamdas.first; lamda2=lamdas.second;
		l_rand=urdist2(rng);
		std::cout<<l_rand<<std::endl;
		lamda=lamda1+(lamda2-lamda1)*l_rand;
		p=get_NewP2(p, lamda, k);
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
	}
	randPoints.push_back(p);
	std::cout<<"full rand point to begin HnR in P: "<<p<<std::endl;
	//sample rnum in Convex body 
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		std::cout<<"rnum is: "<<rnum<<std::endl;
	}
	std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
	
	Point_d p_temp;
	double current_dist, max_dist=0;
	for(std::vector<Point_d>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
		p_temp=Subtrack_points2(*pit,c);
		current_dist=metro_sq(p_temp);
		//current_dist=std::sqrt(current_dist);
		//current_dist=(*pit-c).squared_length();
		if(current_dist>max_dist){
			max_dist=current_dist;
		}
	}
	max_dist=std::sqrt(max_dist);
	std::cout<<"max dist is: "<<max_dist<<std::endl;
	
	int nb1 = dim * (std::log(radius)/std::log(2.0));
    int nb2 = std::ceil(dim * (std::log(max_dist)/std::log(2.0)));
	std::cout<<"nb1: "<<nb1<<"nb2: "<<nb2<<std::endl;
	
	std::vector<double> all_radius;
	
	for(int i=nb1; i<=nb2; ++i){
        all_radius.push_back(std::pow(2.0,(double)i/(double)dim));
        //if (print) {
        //    std::vector<Ball>::iterator bit=balls.end();--bit;
        //    std::cout<<"ball "<<bit-balls.begin()<<" | "<<i
        //           <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
        //}
    }
	
	EXACT_NT telescopic_prod = EXACT_NT(1.0);
	int num_balls=all_radius.size();
	//int nump_PBSmall=10000;
	const double pi = boost::math::constants::pi<double>();
	//double volinit=(2*std::pow(pi,dim/2.0)*std::pow(radius,dim))
     //       / (std::tgamma(dim/2.0)*dim);
       //     radSmall=1.0;
            double vol=0.0;
	for(i=num_balls-1; i>0; i--){
		//if(((double)(rnum))/((double)(nump_PBSmall))>2){
		//	all_radius[0]=radSmall;
		//	break;
		//}
			//const double pi = boost::math::constants::pi<double>();
			vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[i],dim))
            / (std::tgamma(dim/2.0)*dim) * CGAL::to_double(telescopic_prod);
            
        std::cout<<"i is: "<<i<<" vol is: "<<vol<<std::endl;
		if (ball_inside){
			radius=all_radius[i];
			break;
			//all_radius[0]=all_radius[i];
			
		}
		if (randPoints.size()<rnum/1.5){
			ball_inside=true;
		}
		ray_onball=true;
			
			
	
			//return vol;
		//}
		
		// choose a point in PBLarge to be used to generate more rand points
        Point_d p_gen = *randPoints.begin();
        radLarge=all_radius[i];
        std::cout<<"Large rad is: "<<radLarge<<std::endl;
        //i--;
        radSmall=all_radius[i-1];
        std::cout<<"Small rad is: "<<radSmall<<std::endl;

        // num of points in PBSmall and PBLarge
        int nump_PBSmall = 0;
        int nump_PBLarge = randPoints.size();
        std::cout<<"number in large is: "<<nump_PBLarge<<std::endl;
        
        //keep the points in randPoints that fall in PBSmall
        std::vector<Point_d>::iterator rpit=randPoints.begin();
        while(rpit!=randPoints.end()){
            if (isin_ball2(c, radSmall,*rpit) == false){//not in
                rpit=randPoints.erase(rpit);
            } else {
                ++nump_PBSmall;
                ++rpit;
            }
        }
        std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
      //nump_PBLarge = randPoints.size();
        for(int i=1; i<=rnum - nump_PBLarge; ++i){
			for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p_gen, c, radLarge, k, true, ray_onball);
			
			ball_inside=(ball_inside && ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p_gen=get_NewP2(p_gen, lamda, k);
			}
			//std::cout<<"belong to ellips: "<<G.IsIn(p_gen)<<std::endl;
			//std::cout<<"belong to sphere: "<<isin_ball(c, radLarge, p_gen)<<std::endl;
			// count and store in randPoints the points fall in PBSmall
            if (isin_ball2(c, radSmall, p_gen) == true){//is in --- needs creating
                randPoints.push_back(p_gen);
                ++nump_PBSmall;
            }
		}
		
		std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
		std::cout<<"ray_onball is: "<<ray_onball<<std::endl;
		std::cout<<"ball_inside is: "<<ball_inside<<std::endl;
		telescopic_prod *= EXACT_NT(rnum)/EXACT_NT(nump_PBSmall);
		std::cout<<"telescopic prod is: "<<CGAL::to_double(telescopic_prod)<<std::endl;
	}
	vol=0.0;
	//const double pi = boost::math::constants::pi<double>();mpfr_t result,pow,base,exp;
         mpfr_t result,pow,base,exp;
        mpfr_init(result);
        mpfr_init(pow);
        mpfr_init(base);
        mpfr_init(exp);
        mpfr_set_ld(result,2.0,GMP_RNDN);

        mpfr_set_ld(base,pi,GMP_RNDN);
        mpfr_set_ld(exp,dim/2.0,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);

        mpfr_set_ld(base,radius,GMP_RNDN);
        mpfr_set_ld(exp,dim,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);
        mpfr_div_d(result,result,std::tgamma(dim/2.0)*dim,GMP_RNDN);
        mpfr_mul_d(result,result,CGAL::to_double(telescopic_prod),GMP_RNDN);

  //  vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[0],dim))
    //        / (std::tgamma(dim/2.0)*dim)
            //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
    //        * CGAL::to_double(telescopic_prod);
	vol=mpfr_get_d(result,GMP_RNDN);
	
	return vol;
}




double VolEsti_ellips4(ellipsoids G, int walk_len, int rnum, std::vector<double> cntr, double RR){
	
	double radius,lamda1,lamda2,lamda,l_rand,radLarge,radSmall;
	int dim=G.getDim(),i,j,k;
	Point_d c;       //center
	std::pair<double,double> lamdas;
	std::pair<Point_d,double> inscribed_ball;
	bool ball_inside=false, ray_onball=true;
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// the random engine with this seed
	RNGType rng(seed);
	//RNGType rng;
	boost::random::uniform_int_distribution<> uidist1(0,dim-1);
	boost::random::uniform_real_distribution<> urdist2(0,1); 
	
	//get the inscribed ball
	//inscribed_ball=rand_inscribed_ball(G);
	//c=inscribed_ball.first;
	//radius=inscribed_ball.second;
	c=Point_d(dim,cntr.begin(),cntr.end());
	radius=RR;
	std::cout<<"center is: "<<c<<std::endl;
	
	//radius=*rad;
	std::cout<<"radius is: "<<radius<<std::endl;
	
	CGAL::Random_points_in_ball_d<Point_d> gen (dim, radius);
    Point_d p = *gen;
    //p = p + (c-CGAL::Origin());
    p=Add_points2(p,c);
    std::vector<Point_d> randPoints;
    std::cout<<"first rand point to begin HnR in P: "<<p<<std::endl;
    std::cout<<G.IsIn(p)<<std::endl;
	
	k = uidist1(rng); //rand coord
	
	//for(i=0; i<
	lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball); //needs changes
	lamda1=lamdas.first;
	lamda2=lamdas.second;
	l_rand=urdist2(rng);
	lamda=lamda1+(lamda2-lamda1)*l_rand;
	p=get_NewP2(p, lamda, k);
	std::cout<<"second rand point to begin HnR in P: "<<p<<std::endl;
	for (i=0; i<1000; i++){
		k = uidist1(rng);
		std::cout<<k<<std::endl;
		lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
		lamda1=lamdas.first; lamda2=lamdas.second;
		l_rand=urdist2(rng);
		std::cout<<l_rand<<std::endl;
		lamda=lamda1+(lamda2-lamda1)*l_rand;
		p=get_NewP2(p, lamda, k);
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
	}
	randPoints.push_back(p);
	std::cout<<"full rand point to begin HnR in P: "<<p<<std::endl;
	//sample rnum in Convex body 
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		std::cout<<"rnum is: "<<rnum<<std::endl;
	}
	std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
	
	
	
	Point_d p_temp;
	double current_dist, max_dist=0;
	for(std::vector<Point_d>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
		p_temp=Subtrack_points2(*pit,c);
		current_dist=metro_sq(p_temp);
		//current_dist=std::sqrt(current_dist);
		//current_dist=(*pit-c).squared_length();
		if(current_dist>max_dist){
			max_dist=current_dist;
		}
	}
	max_dist=std::sqrt(max_dist);
	std::cout<<"max dist is: "<<max_dist<<std::endl;
	/*randPoints.clear();
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		//std::cout<<"rnum is: "<<rnum<<std::endl;
	}*/
	
	int nb1 = dim * (std::log(radius)/std::log(2.0));
    int nb2 = std::ceil(dim * (std::log(max_dist)/std::log(2.0)));
	std::cout<<"nb1: "<<nb1<<"nb2: "<<nb2<<std::endl;
	
	std::vector<double> all_radius;
	
	for(int i=nb1; i<=nb2; ++i){
        all_radius.push_back(std::pow(2.0,(double)i/(double)dim));
        //if (print) {
        //    std::vector<Ball>::iterator bit=balls.end();--bit;
        //    std::cout<<"ball "<<bit-balls.begin()<<" | "<<i
        //           <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
        //}
    }
	
	EXACT_NT telescopic_prod = EXACT_NT(1.0);
	int num_balls=all_radius.size();
	//int nump_PBSmall=10000;
	const double pi = boost::math::constants::pi<double>();
	//double volinit=(2*std::pow(pi,dim/2.0)*std::pow(radius,dim))
     //       / (std::tgamma(dim/2.0)*dim);
       //     radSmall=1.0;
            double vol=0.0;
	for(i=num_balls-1; i>0; i--){
		//if(((double)(rnum))/((double)(nump_PBSmall))>2){
		//	all_radius[0]=radSmall;
		//	break;
		//}
			//const double pi = boost::math::constants::pi<double>();
			vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[i],dim))
            / (std::tgamma(dim/2.0)*dim) * CGAL::to_double(telescopic_prod);
            
        std::cout<<"i is: "<<i<<" vol is: "<<vol<<std::endl;
		if (ball_inside){
			
			radius=all_radius[i];
			break;
			//all_radius[0]=all_radius[i];
			
		}
		if (randPoints.size()<rnum/1.5){
			ball_inside=true;
		}
		ray_onball=true;
			
			
	
			//return vol;
		//}
		
		// choose a point in PBLarge to be used to generate more rand points
        Point_d p_gen = *randPoints.begin();
        radLarge=all_radius[i];
        std::cout<<"Large rad is: "<<radLarge<<std::endl;
        //i--;
        radSmall=all_radius[i-1];
        std::cout<<"Small rad is: "<<radSmall<<std::endl;

        // num of points in PBSmall and PBLarge
        int nump_PBSmall = 0;
        int nump_PBLarge = randPoints.size();
        std::cout<<"number in large is: "<<nump_PBLarge<<std::endl;
        
        //keep the points in randPoints that fall in PBSmall
        std::vector<Point_d>::iterator rpit=randPoints.begin();
        while(rpit!=randPoints.end()){
            if (isin_ball2(c, radSmall,*rpit) == false){//not in
                rpit=randPoints.erase(rpit);
            } else {
                ++nump_PBSmall;
                ++rpit;
            }
        }
        std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
      //nump_PBLarge = randPoints.size();
        for(int i=1; i<=rnum - nump_PBLarge; ++i){
			for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p_gen, c, radLarge, k, true, ray_onball);
			
			ball_inside=(ball_inside && ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p_gen=get_NewP2(p_gen, lamda, k);
			}
			//std::cout<<"belong to ellips: "<<G.IsIn(p_gen)<<std::endl;
			//std::cout<<"belong to sphere: "<<isin_ball(c, radLarge, p_gen)<<std::endl;
			// count and store in randPoints the points fall in PBSmall
            if (isin_ball2(c, radSmall, p_gen) == true){//is in --- needs creating
                randPoints.push_back(p_gen);
                ++nump_PBSmall;
            }
		}
		
		std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
		std::cout<<"ray_onball is: "<<ray_onball<<std::endl;
		std::cout<<"ball_inside is: "<<ball_inside<<std::endl;
		telescopic_prod *= EXACT_NT(rnum)/EXACT_NT(nump_PBSmall);
		std::cout<<"telescopic prod is: "<<CGAL::to_double(telescopic_prod)<<std::endl;
	}
	vol=0.0;
	//const double pi = boost::math::constants::pi<double>();mpfr_t result,pow,base,exp;
         mpfr_t result,pow,base,exp;
        mpfr_init(result);
        mpfr_init(pow);
        mpfr_init(base);
        mpfr_init(exp);
        mpfr_set_ld(result,2.0,GMP_RNDN);

        mpfr_set_ld(base,pi,GMP_RNDN);
        mpfr_set_ld(exp,dim/2.0,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);

        mpfr_set_ld(base,radius,GMP_RNDN);
        mpfr_set_ld(exp,dim,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);
        mpfr_div_d(result,result,std::tgamma(dim/2.0)*dim,GMP_RNDN);
        mpfr_mul_d(result,result,CGAL::to_double(telescopic_prod),GMP_RNDN);

  //  vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[0],dim))
    //        / (std::tgamma(dim/2.0)*dim)
            //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
    //        * CGAL::to_double(telescopic_prod);
	vol=mpfr_get_d(result,GMP_RNDN);
	
	return vol;
}


//volesti: unit simplex, ellipsoid and two parallel hyperplanes
double VolEsti_ellips5(ellipsoids G, int walk_len, int rnum, std::vector<double> facet, double z1, double z2){
	
	double radius,lamda1,lamda2,lamda,l_rand,radLarge,radSmall;
	int dim=G.getDim(),i,j,k;
	Point_d c;       //center
	std::pair<double,double> lamdas;
	std::pair<Point_d,double> inscribed_ball;
	bool ball_inside=false, ray_onball=true;
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// the random engine with this seed
	RNGType rng(seed);
	//RNGType rng;
	boost::random::uniform_int_distribution<> uidist1(0,dim-1);
	boost::random::uniform_real_distribution<> urdist2(0,1); 
	
	//get the inscribed ball
	inscribed_ball=rand_inscribed_ball2(G,facet,z1,z2);
	c=inscribed_ball.first;
	radius=inscribed_ball.second;
	//c=Point_d(dim,cntr.begin(),cntr.end());
	//radius=RR;
	std::cout<<"center is: "<<c<<std::endl;
	
	//radius=*rad;
	std::cout<<"radius is: "<<radius<<std::endl;
	
	CGAL::Random_points_in_ball_d<Point_d> gen (dim, radius);
    Point_d p = *gen;
    //p = p + (c-CGAL::Origin());
    p=Add_points2(p,c);
    std::vector<Point_d> randPoints;
    std::cout<<"first rand point to begin HnR in P: "<<p<<std::endl;
    std::cout<<G.IsIn(p)<<std::endl;
	
	k = uidist1(rng); //rand coord
	
	//for(i=0; i<
	lamdas=hit_and_run_newP3(G, p, c, radius, facet, z1, z2, k, false, ray_onball); //needs changes
	lamda1=lamdas.first;
	lamda2=lamdas.second;
	l_rand=urdist2(rng);
	lamda=lamda1+(lamda2-lamda1)*l_rand;
	p=get_NewP2(p, lamda, k);
	std::cout<<"second rand point to begin HnR in P: "<<p<<std::endl;
	for (i=0; i<1000; i++){
		k = uidist1(rng);
		std::cout<<k<<std::endl;
		lamdas=hit_and_run_newP3(G, p, c, radius, facet, z1, z2, k, false, ray_onball);
		lamda1=lamdas.first; lamda2=lamdas.second;
		l_rand=urdist2(rng);
		std::cout<<l_rand<<std::endl;
		lamda=lamda1+(lamda2-lamda1)*l_rand;
		p=get_NewP2(p, lamda, k);
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
	}
	randPoints.push_back(p);
	std::cout<<"full rand point to begin HnR in P: "<<p<<std::endl;
	//sample rnum in Convex body 
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP3(G, p, c, radius, facet, z1, z2, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		//std::cout<<"rnum is: "<<rnum<<std::endl;
	}
	std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
	
	
	
	Point_d p_temp;
	double current_dist, max_dist=0;
	for(std::vector<Point_d>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
		p_temp=Subtrack_points2(*pit,c);
		current_dist=metro_sq(p_temp);
		//current_dist=std::sqrt(current_dist);
		//current_dist=(*pit-c).squared_length();
		if(current_dist>max_dist){
			max_dist=current_dist;
		}
	}
	max_dist=std::sqrt(max_dist);
	std::cout<<"max dist is: "<<max_dist<<std::endl;
	/*randPoints.clear();
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		//std::cout<<"rnum is: "<<rnum<<std::endl;
	}*/
	
	int nb1 = dim * (std::log(radius)/std::log(2.0));
    int nb2 = std::ceil(dim * (std::log(max_dist)/std::log(2.0)));
	std::cout<<"nb1: "<<nb1<<"nb2: "<<nb2<<std::endl;
	
	std::vector<double> all_radius;
	
	for(int i=nb1; i<=nb2; ++i){
        all_radius.push_back(std::pow(2.0,(double)i/(double)dim));
        //if (print) {
        //    std::vector<Ball>::iterator bit=balls.end();--bit;
        //    std::cout<<"ball "<<bit-balls.begin()<<" | "<<i
        //           <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
        //}
    }
	
	EXACT_NT telescopic_prod = EXACT_NT(1.0);
	int num_balls=all_radius.size();
	//int nump_PBSmall=10000;
	const double pi = boost::math::constants::pi<double>();
	//double volinit=(2*std::pow(pi,dim/2.0)*std::pow(radius,dim))
     //       / (std::tgamma(dim/2.0)*dim);
       //     radSmall=1.0;
            double vol=0.0;
	for(i=num_balls-1; i>0; i--){
		//if(((double)(rnum))/((double)(nump_PBSmall))>2){
		//	all_radius[0]=radSmall;
		//	break;
		//}
			//const double pi = boost::math::constants::pi<double>();
			vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[i],dim))
            / (std::tgamma(dim/2.0)*dim) * CGAL::to_double(telescopic_prod);
            
        std::cout<<"i is: "<<i<<" vol is: "<<vol<<std::endl;
		if (ball_inside){
			
			radius=all_radius[i];
			break;
			//all_radius[0]=all_radius[i];
			
		}
		if (randPoints.size()<rnum/1.5){
			ball_inside=true;
		}
		ray_onball=true;
			
			
	
			//return vol;
		//}
		
		// choose a point in PBLarge to be used to generate more rand points
        Point_d p_gen = *randPoints.begin();
        radLarge=all_radius[i];
        std::cout<<"Large rad is: "<<radLarge<<std::endl;
        //i--;
        radSmall=all_radius[i-1];
        std::cout<<"Small rad is: "<<radSmall<<std::endl;

        // num of points in PBSmall and PBLarge
        int nump_PBSmall = 0;
        int nump_PBLarge = randPoints.size();
        std::cout<<"number in large is: "<<nump_PBLarge<<std::endl;
        
        //keep the points in randPoints that fall in PBSmall
        std::vector<Point_d>::iterator rpit=randPoints.begin();
        while(rpit!=randPoints.end()){
            if (isin_ball2(c, radSmall,*rpit) == false){//not in
                rpit=randPoints.erase(rpit);
            } else {
                ++nump_PBSmall;
                ++rpit;
            }
        }
        std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
      //nump_PBLarge = randPoints.size();
        for(int i=1; i<=rnum - nump_PBLarge; ++i){
			for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP3(G, p_gen, c, radLarge, facet, z1, z2, k, true, ray_onball);
			
			ball_inside=(ball_inside && ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p_gen=get_NewP2(p_gen, lamda, k);
			}
			//std::cout<<"belong to ellips: "<<G.IsIn(p_gen)<<std::endl;
			//std::cout<<"belong to sphere: "<<isin_ball(c, radLarge, p_gen)<<std::endl;
			// count and store in randPoints the points fall in PBSmall
            if (isin_ball2(c, radSmall, p_gen) == true){//is in --- needs creating
                randPoints.push_back(p_gen);
                ++nump_PBSmall;
            }
		}
		
		std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
		std::cout<<"ray_onball is: "<<ray_onball<<std::endl;
		std::cout<<"ball_inside is: "<<ball_inside<<std::endl;
		telescopic_prod *= EXACT_NT(rnum)/EXACT_NT(nump_PBSmall);
		std::cout<<"telescopic prod is: "<<CGAL::to_double(telescopic_prod)<<std::endl;
	}
	vol=0.0;
	//const double pi = boost::math::constants::pi<double>();mpfr_t result,pow,base,exp;
         mpfr_t result,pow,base,exp;
        mpfr_init(result);
        mpfr_init(pow);
        mpfr_init(base);
        mpfr_init(exp);
        mpfr_set_ld(result,2.0,GMP_RNDN);

        mpfr_set_ld(base,pi,GMP_RNDN);
        mpfr_set_ld(exp,dim/2.0,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);

        mpfr_set_ld(base,radius,GMP_RNDN);
        mpfr_set_ld(exp,dim,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);
        mpfr_div_d(result,result,std::tgamma(dim/2.0)*dim,GMP_RNDN);
        mpfr_mul_d(result,result,CGAL::to_double(telescopic_prod),GMP_RNDN);

  //  vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[0],dim))
    //        / (std::tgamma(dim/2.0)*dim)
            //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
    //        * CGAL::to_double(telescopic_prod);
	vol=mpfr_get_d(result,GMP_RNDN);
	
	return vol;
}


double VolEsti_ellips_nonConvex(ellipsoids G1, ellipsoids G2, int walk_len, int rnum, std::vector<double> facet, double z1, double z2){
	
	double radius,lamda1,lamda2,lamda,l_rand,radLarge,radSmall;
	int dim=G1.getDim(),i,j,k;
	Point_d c;       //center
	std::pair<double,double> lamdas;
	std::pair<Point_d,double> inscribed_ball;
	bool ball_inside=false, ray_onball=true;
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// the random engine with this seed
	RNGType rng(seed);
	//RNGType rng;
	boost::random::uniform_int_distribution<> uidist1(0,dim-1);
	boost::random::uniform_real_distribution<> urdist2(0,1); 
	
	//get the inscribed ball
	inscribed_ball=rand_inscribed_ball_nonConv(G1,G2,facet,z1,z2);
	c=inscribed_ball.first;
	radius=inscribed_ball.second;
	//c=Point_d(dim,cntr.begin(),cntr.end());
	//radius=RR;
	std::cout<<"center is: "<<c<<std::endl;
	
	//radius=*rad;
	std::cout<<"radius is: "<<radius<<std::endl;
	
	CGAL::Random_points_in_ball_d<Point_d> gen (dim, radius);
    Point_d p = *gen;
    //p = p + (c-CGAL::Origin());
    p=Add_points2(p,c);
    std::vector<Point_d> randPoints;
    std::cout<<"first rand point to begin HnR in P: "<<p<<std::endl;
    std::cout<<G1.IsIn(p)<<std::endl;
	
	k = uidist1(rng); //rand coord
	
	//for(i=0; i<
	lamdas=hit_and_run_newP_nonConv(G1, G2, p, c, radius, facet, z1, z2, k, false, ray_onball); //needs changes
	lamda1=lamdas.first;
	lamda2=lamdas.second;
	l_rand=urdist2(rng);
	lamda=lamda1+(lamda2-lamda1)*l_rand;
	p=get_NewP2(p, lamda, k);
	std::cout<<"second rand point to begin HnR in P: "<<p<<std::endl;
	for (i=0; i<1000; i++){
		k = uidist1(rng);
		std::cout<<k<<std::endl;
		lamdas=hit_and_run_newP_nonConv(G1, G2, p, c, radius, facet, z1, z2, k, false, ray_onball);
		lamda1=lamdas.first; lamda2=lamdas.second;
		l_rand=urdist2(rng);
		std::cout<<l_rand<<std::endl;
		lamda=lamda1+(lamda2-lamda1)*l_rand;
		p=get_NewP2(p, lamda, k);
		std::cout<<G1.IsIn(p)<<" "<<G2.IsIn(p)<<" first"<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
	}
	randPoints.push_back(p);
	std::cout<<"full rand point to begin HnR in P: "<<p<<std::endl;
	//sample rnum in Convex body 
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP_nonConv(G1, G2, p, c, radius, facet, z1, z2, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		std::cout<<"rnum is: "<<rnum<<std::endl;
		//std::cout<<G1.IsIn(c)<<" "<<G2.IsIn(c)<<std::endl;
	}
	std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
	
	
	
	Point_d p_temp;
	double current_dist, max_dist=0;
	for(std::vector<Point_d>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
		p_temp=Subtrack_points2(*pit,c);
		current_dist=metro_sq(p_temp);
		//current_dist=std::sqrt(current_dist);
		//current_dist=(*pit-c).squared_length();
		if(current_dist>max_dist){
			max_dist=current_dist;
		}
	}
	max_dist=std::sqrt(max_dist);
	std::cout<<"max dist is: "<<max_dist<<std::endl;
	/*randPoints.clear();
	for (i=0; i<rnum; i++){
		for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			lamdas=hit_and_run_newP2(G, p, c, radius, k, false, ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p=get_NewP2(p, lamda, k);
		}
		randPoints.push_back(p);
		//std::cout<<"number in rand points in P: "<<randPoints.size()<<std::endl;
		//std::cout<<"belong to ellips: "<<G.IsIn(p)<<std::endl;
		//std::cout<<"belong to sphere: "<<isin_ball(c, radius, p)<<std::endl;
		//std::cout<<"rnum is: "<<rnum<<std::endl;
	}*/
	
	int nb1 = dim * (std::log(radius)/std::log(2.0));
    int nb2 = std::ceil(dim * (std::log(max_dist)/std::log(2.0)));
	std::cout<<"nb1: "<<nb1<<"nb2: "<<nb2<<std::endl;
	
	std::vector<double> all_radius;
	
	for(int i=nb1; i<=nb2; ++i){
        all_radius.push_back(std::pow(2.0,(double)i/(double)dim));
        //if (print) {
        //    std::vector<Ball>::iterator bit=balls.end();--bit;
        //    std::cout<<"ball "<<bit-balls.begin()<<" | "<<i
        //           <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
        //}
    }
	
	EXACT_NT telescopic_prod = EXACT_NT(1.0);
	int num_balls=all_radius.size();
	//int nump_PBSmall=10000;
	const double pi = boost::math::constants::pi<double>();
	//double volinit=(2*std::pow(pi,dim/2.0)*std::pow(radius,dim))
     //       / (std::tgamma(dim/2.0)*dim);
       //     radSmall=1.0;
            double vol=0.0;
	for(i=num_balls-1; i>0; i--){
		//if(((double)(rnum))/((double)(nump_PBSmall))>2){
		//	all_radius[0]=radSmall;
		//	break;
		//}
			//const double pi = boost::math::constants::pi<double>();
			vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[i],dim))
            / (std::tgamma(dim/2.0)*dim) * CGAL::to_double(telescopic_prod);
            
        std::cout<<"i is: "<<i<<" vol is: "<<vol<<std::endl;
		if (ball_inside){
			
			radius=all_radius[i];
			break;
			//all_radius[0]=all_radius[i];
			
		}
		if (randPoints.size()<rnum){
			ball_inside=true;
		}
		ray_onball=true;
			
			
	
			//return vol;
		//}
		
		// choose a point in PBLarge to be used to generate more rand points
        Point_d p_gen = *randPoints.begin();
        radLarge=all_radius[i];
        std::cout<<"Large rad is: "<<radLarge<<std::endl;
        //i--;
        radSmall=all_radius[i-1];
        std::cout<<"Small rad is: "<<radSmall<<std::endl;

        // num of points in PBSmall and PBLarge
        int nump_PBSmall = 0;
        int nump_PBLarge = randPoints.size();
        std::cout<<"number in large is: "<<nump_PBLarge<<std::endl;
        
        //keep the points in randPoints that fall in PBSmall
        std::vector<Point_d>::iterator rpit=randPoints.begin();
        while(rpit!=randPoints.end()){
            if (isin_ball2(c, radSmall,*rpit) == false){//not in
                rpit=randPoints.erase(rpit);
            } else {
                ++nump_PBSmall;
                ++rpit;
            }
        }
        std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
      //nump_PBLarge = randPoints.size();
        for(int i=1; i<=rnum - nump_PBLarge; ++i){
			for (j=0; j<walk_len; j++){
			k = uidist1(rng);
			//lamdas=hit_and_run_newP3(G, p_gen, c, radLarge, facet, z1, z2, k, true, ray_onball);
			lamdas=hit_and_run_newP_nonConv(G1, G2, p_gen, c,radLarge, facet, z1, z2, k, true, ray_onball);
			
			ball_inside=(ball_inside && ray_onball);
			lamda1=lamdas.first; lamda2=lamdas.second;
			l_rand=urdist2(rng);
			lamda=lamda1+(lamda2-lamda1)*l_rand;
			p_gen=get_NewP2(p_gen, lamda, k);
			}
			//std::cout<<"belong to ellips: "<<G.IsIn(p_gen)<<std::endl;
			//std::cout<<"belong to sphere: "<<isin_ball(c, radLarge, p_gen)<<std::endl;
			// count and store in randPoints the points fall in PBSmall
            if (isin_ball2(c, radSmall, p_gen) == true){//is in --- needs creating
                randPoints.push_back(p_gen);
                ++nump_PBSmall;
            }
		}
		
		std::cout<<"number in small is: "<<nump_PBSmall<<std::endl;
		std::cout<<"ray_onball is: "<<ray_onball<<std::endl;
		std::cout<<"ball_inside is: "<<ball_inside<<std::endl;
		telescopic_prod *= EXACT_NT(rnum)/EXACT_NT(nump_PBSmall);
		//std::cout<<"telescopic prod is: "<<CGAL::to_double(telescopic_prod)<<std::endl;
	}
	vol=0.0;
	//const double pi = boost::math::constants::pi<double>();mpfr_t result,pow,base,exp;
         mpfr_t result,pow,base,exp;
        mpfr_init(result);
        mpfr_init(pow);
        mpfr_init(base);
        mpfr_init(exp);
        mpfr_set_ld(result,2.0,GMP_RNDN);

        mpfr_set_ld(base,pi,GMP_RNDN);
        mpfr_set_ld(exp,dim/2.0,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);

        mpfr_set_ld(base,radius,GMP_RNDN);
        mpfr_set_ld(exp,dim,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);
        mpfr_div_d(result,result,std::tgamma(dim/2.0)*dim,GMP_RNDN);
        mpfr_mul_d(result,result,CGAL::to_double(telescopic_prod),GMP_RNDN);

  //  vol = (2*std::pow(pi,dim/2.0)*std::pow(all_radius[0],dim))
    //        / (std::tgamma(dim/2.0)*dim)
            //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
    //        * CGAL::to_double(telescopic_prod);
	vol=mpfr_get_d(result,GMP_RNDN);
	
	return vol;
}



//----------------------


double sample_2hyp_par(int dim, int num,std::vector<double> pl, double z1, double z2){
	
	std::vector<Point_d> points;
	Point_d p;
	int sum=0,i,j;
	double sum_p;
	
	Sam_Unit(dim, num, points);
	
	for (i=0; i<num; i++){
		p=points[i];
		sum_p=0.0;
		for (j=0; j<dim; j++){
			sum_p+=p[j]*pl[j];
		}
		if (sum_p<z2 && sum_p>z1){
			sum++;
		}
	}
	std::cout<<"sum is: "<<sum<<std::endl;
	return ((double)sum)/((double)num);
}


double sample_2hyp(int dim, int num, Plane_d pl1, Plane_d pl2){
	
	std::vector<Point_d> points;
	Point_d p;
	int sum=0,i,j;
	double sum_p1,sum_p2;
	
	Sam_Unit(dim, num, points);
	
	for (i=0; i<num; i++){
		p=points[i];
		sum_p1=0.0, sum_p2=0.0;
		for (j=0; j<dim; j++){
			sum_p1+=p[j]*pl1[j];
            sum_p2+=p[j]*pl2[j];
            
        }
		if (sum_p2<(pl2[dim]) && sum_p1>(pl1[dim])){
			sum++;
		}
	}
	std::cout<<"sum is: "<<sum<<std::endl;
	return ((double)sum)/((double)num);
}


double sample_4hyp_par(int dim, int num,std::vector<double> pl1, double z11, double z12, std::vector<double> pl2, double z21, double z22){
	
	std::vector<Point_d> points;
	Point_d p;
	int sum=0,i,j;
	double sum_p1,sum_p2;
	
	Sam_Unit(dim, num, points);
	
	for (i=0; i<num; i++){
		p=points[i];
		sum_p1=0.0; sum_p2=0.0;
		for (j=0; j<dim; j++){
			sum_p1+=p[j]*pl1[j];
			sum_p2+=p[j]*pl2[j];
		}
		if (sum_p1<z12 && sum_p1>z11 && sum_p2<z22 && sum_p2>z21){
			sum++;
		}
	}
	std::cout<<"sum is: "<<sum<<std::endl;
	return ((double)sum)/((double)num);
}



double sample_cut_ellipsoid(int dim, int num, ellipsoids G){
	
	std::vector<Point_d> points;
	Point_d p;
	int sum=0,i;
	
	Sam_Unit(dim, num, points);
	
	for (i=0; i<num; i++){
		p=points[i];
		if (G.IsIn(p)){
			sum++;
		}
	}
	std::cout<<"sum is: "<<sum<<std::endl;
	return (double)(sum)/(double)(num);
	
	
	
	
	
}



double sample_cut_ellipsoid_and_hyps(int dim, int num, ellipsoids G, std::vector<double> pl, double z1, double z2){
	
	std::vector<Point_d> points;
	Point_d p;
	int sum=0,i,j;
	double sum_p;
	
	Sam_Unit(dim, num, points);
	
	for (i=0; i<num; i++){
		p=points[i]; sum_p=0.0;
		if (G.IsIn(p)){
			for (j=0; j<dim; j++){
				sum_p+=p[j]*pl[j];
			}
			if (sum_p<z2 && sum_p>z1){
				sum++;
			}
		}
	}
	std::cout<<"sum is: "<<sum<<std::endl;
	return (double)(sum)/(double)(num);
	
	
	
	
	
}



double sampleToArb(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end,Plane_d plane1, Plane_d plane2, int num){

	int i,j,sum=0,dim=std::distance(it_beg,it_end);
	double side1,side2,vol;
	Point_d p, p0=*it_beg;;
	std::vector<Point_d> points;
	std::vector<Point_d>::iterator iter;
	dim--;
	
	MatrixXd A(dim,dim);
	
	for (j=1; j<dim+1; j++){
		Point_d pk=*(it_beg+j);
		for (i=0; i<dim; i++){			
			A(i,j-1)=pk[i]-p0[i];
		}
	}
	vol=A.determinant();
	
	Sam_arbest(it_beg,it_end,num,points);
	//sam_simplex2(it_beg,it_end,num,points);

	for (iter=points.begin(); iter!=points.end(); iter++){
		p=*iter;
		side1=0.0; side2=0.0;
		for (j=0; j<dim; j++){
			side1+=plane1[j]*p[j];
			side2+=plane2[j]*p[j];
		}
		//side+=plane[dim];
		if (side1<(-plane1[dim]) && side2<(-plane2[dim])){
			sum+=1;
		}
	}
	//std::cout<<"Number of Points in Arbitrary Simplex: "<<sum<<std::endl;
	vol= ((double)sum)/((double)num);
	for (i=0; i<dim; i++){
		vol=vol/((double)(i+1));
	}
	
	return vol;

}


double sample_cut_ellipsoid_and_hyps_nonConv(int dim, int num, ellipsoids G1, ellipsoids G2, std::vector<double> pl, double z1, double z2){
	
	std::vector<Point_d> points;
	Point_d p;
	int sum=0,i,j;
	double sum_p;
	
	Sam_Unit(dim, num, points);
	
	for (i=0; i<num; i++){
		p=points[i]; sum_p=0.0;
		if (G1.IsIn(p) && !G2.IsIn(p)){
			for (j=0; j<dim; j++){
				sum_p+=p[j]*pl[j];
			}
			if (sum_p<z2 && sum_p>z1){
				sum++;
			}
		}
	}
	std::cout<<"sum is: "<<sum<<std::endl;
	return (double)(sum)/(double)(num);
	
	
	
	
	
}


double sampleToUnit(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end,Plane_d plane1,Plane_d plane2, int num){
	

	int j,i,sum;
	Point_d p0=*it_beg,p1;
	int dim=p0.dimension();
	double side1,side2,volume,z1,z2;
	bool t=false;
	//int pivot11[dim],pivot12[dim],pivot21[dim],pivot22[dim];
	std::vector<Point_d> points;
	std::vector<Point_d>::iterator iter1;
	
	MatrixXd A(dim,dim);
	VectorXd pl1(dim);
	VectorXd pl2(dim);
	
	for (i=0; i<dim; i++){
		pl1(i)=plane1[i];
		pl2(i)=plane2[i];
	}
	z1=-plane1[dim];
	z2=-plane2[dim];
	
	for (j=1; j<dim+1; j++){
		Point_d pk=*(it_beg+j);
		z1=z1-pl1(j-1)*p0[j-1];
		z2=z2-pl2(j-1)*p0[j-1];
		for (i=0; i<dim; i++){			
			A(i,j-1)=pk[i]-p0[i];
		}
	}
	pl1=pl1.transpose()*A;
	pl2=pl2.transpose()*A;
	
	Sam_Unit(dim,num,points);
	volume=A.determinant();
	sum=0;
	 
	for (i=0; i<dim; i++){
		p1=points[i];
		side1=0.0; side2=0.0;
		for (j=0; j<dim; j++){
			side1+=pl1(j)*p1[j];
			side2+=pl2(j)*p1[j];
		}
		
		if (side1<z1 && side2<z2){
			sum++;
		}
	}
	volume=volume*( ((double)sum)/((double)num) );
	
	for (i=0; i<dim; i++){
		volume=volume/((double)(i+1));
	}
	
	return volume;
	
}
	

std::vector<std::vector<double> > get_par_4hyp_vols(int dim, int num, std::vector<double>  pl1, std::vector<double> pl2){
	
	int i,j,col,row;
	std::vector<double> vec1,vec2,Zs1,Zs2;
	double sum1,sum2,sum;
	std::vector<Point_d> points;
	std::pair< std::vector<double>,std::vector<double> > result;
	Point_d p;
	
	Sam_Unit_NoProjection(dim, num, points);
	
    std::vector<std::vector<int> > Matrix(100);
    std::vector<std::vector<double> > pos_Matrix(100);
	for (i=0; i<100; i++){
		Matrix[i].resize(100);
        pos_Matrix[i].resize(100);
	}
	for (i=0; i<100; i++){
		for (j=0; j<100; j++){
			Matrix[i][j]=0;
		}
	}
    
	for (i=0; i<num; i++){
		p=points[i];
		//std::cout<<p<<std::endl;
		sum1=0.0; sum2=0.0; sum=0.0;
		for (j=0; j<dim; j++){
			sum1+=p[j]*pl1[j];
			sum2+=p[j]*pl2[j];
			//sum+=p[j];
			
		}
		//std::cout<<sum<<std::endl;
		vec1.push_back(sum1);
		vec2.push_back(sum2);
	}
	std::sort( vec1.begin(), vec1.end() );
	std::sort( vec2.begin(), vec2.end() );
	
	for (i=1; i<100; i++){
		//std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
		//std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
		Zs1.push_back(vec1[((int)std::floor(i*(0.01)*((double)num) ))]);
		Zs2.push_back(vec2[((int)std::floor(i*(0.01)*((double)num) ))]);
	}
	
	
	
	//std::cout<<"hello2"<<std::endl;
	for (i=0; i<num; i++){
		p=points[i];
		//std::cout<<"dimension is: "<<p.dimension()<<std::endl;
		sum1=0.0; sum2=0.0;
		col=-1; row=-1;
		for (j=0; j<dim; j++){
			sum1+=p[j]*pl1[j];
			sum2+=p[j]*pl2[j];
		}
		//std::cout<<"hello3"<<std::endl;
		for (j=0; j<Zs1.size(); j++){
			if (sum1<Zs1[j]){
				col=j;
				break;
			}
		}
		for (j=0; j<Zs2.size(); j++){
			if (sum2<Zs2[j]){
				row=j;
				break;
			}
		}
		if (col==-1){
			col=99;
		}
		
		if (row==-1){
			row=99;
		}
		//std::cout<<"hello4"<<std::endl;
		Matrix[row][col]++;
	}
	
	for (i=0; i<100; i++){
		for (j=0; j<100; j++){
			pos_Matrix[i][j]=((double)Matrix[i][j])/((double)num);
		}
	}
	
	
	return pos_Matrix;
	
}



std::vector<std::vector<double> > get_par_hypellips_vols(int dim, int num, std::vector<double>  pl, ellipsoids G){
	
	int i,j,col,row;
	std::vector<double> vec1,vec2,Zs1,Cs;
	double sum1,sum2,sum;
	std::vector<Point_d> points;
	std::pair< std::vector<double>,std::vector<double> > result;
	Point_d p;
	
	Sam_Unit_NoProjection(dim, num, points);
	
    std::vector<std::vector<int> > Matrix(100);
    std::vector<std::vector<double> > pos_Matrix(100);
	for (i=0; i<100; i++){
		Matrix[i].resize(100);
        pos_Matrix[i].resize(100);
	}
	for (i=0; i<100; i++){
		for (j=0; j<100; j++){
			Matrix[i][j]=0;
		}
	}
    
	for (i=0; i<num; i++){
		p=points[i];
		//std::cout<<p<<std::endl;
		sum1=0.0;
        sum2=G.MatMult(p);
		for (j=0; j<dim; j++){
			sum1+=p[j]*pl[j];
			//sum2+=p[j]*pl2[j];
			//sum+=p[j];
			
		}
		//std::cout<<sum<<std::endl;
		vec1.push_back(sum1);
		vec2.push_back(sum2);
	}
	std::sort( vec1.begin(), vec1.end() );
	std::sort( vec2.begin(), vec2.end() );
	
	for (i=1; i<100; i++){
		//std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
		//std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
		Zs1.push_back(vec1[((int)std::floor(i*(0.01)*((double)num) ))]);
		Cs.push_back(vec2[((int)std::floor(i*(0.01)*((double)num) ))]);
	}
	
	
	
	//std::cout<<"hello2"<<std::endl;
	for (i=0; i<num; i++){
		p=points[i];
		//std::cout<<"dimension is: "<<p.dimension()<<std::endl;
		sum1=0.0; sum2=0.0;
		col=-1; row=-1;
        sum2=G.MatMult(p);
		for (j=0; j<dim; j++){
			sum1+=p[j]*pl[j];
			//sum2+=p[j]*pl2[j];
		}
		//std::cout<<"hello3"<<std::endl;
		for (j=0; j<Zs1.size(); j++){
			if (sum1<Zs1[j]){
				col=j;
				break;
			}
		}
		for (j=0; j<Cs.size(); j++){
			if (sum2<Cs[j]){
				row=j;
				break;
			}
		}
		if (col==-1){
			col=99;
		}
		
		if (row==-1){
			row=99;
		}
		//std::cout<<"hello4"<<std::endl;
		Matrix[row][col]++;
	}
	
	for (i=0; i<100; i++){
		for (j=0; j<100; j++){
			pos_Matrix[i][j]=((double)Matrix[i][j])/((double)num);
		}
	}
	
	
	return pos_Matrix;
	
}


