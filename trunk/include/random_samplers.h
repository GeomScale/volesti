// RandGeom is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// RandGeom is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.
// 
// Developer: Vissarion Fisikopoulos

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H

// function to find intersection of a line and a polytope 
/*OLD
template <class T>
Vector line_bisect(Point pin, 
                      Vector l, 
                      T &P, 
											vars &var){
	
	double err = var.err;											
  Vector vin=pin-CGAL::Origin();
  //first compute a point outside P along the line
  Point pout=pin;
  //std::cout<<"Starting inside point: "<<vin<<std::endl;
	
	Vector aug(l);
  while(P.is_in(pout) == -1){
    aug*=2;
    pout+=aug;
    //std::cout<<"l="<<l<<" aug="<<aug<<" Outside point:(pout) "<<pout<<std::endl;
  }
  Vector vout=pout-CGAL::Origin();
  
  //intersect using bisection
  //std::cout<<"pout"<<vout<<std::endl;
  Vector vmid;
  double len;
  do{
		vmid=(vin+vout)/2;
		if(P.is_in(CGAL::Origin()+vmid) != -1)
			vout=vmid;
		else
			vin=vmid;
		len=CGAL::to_double((vin-vout).squared_length());
		//std::cout<<"len="<<bool(len<err)<<std::endl;
	}while(len > err);
  
  //std::cout<<"Intersection point: ";
  //round_print(vmid);
  //return vmid; 
	
	return vin; //ensure that the point is always in P 
}
*/

/* Hit and run with random directions */
/*OLD
template <class T>
int hit_and_run(Point &p,
					      T &P,
					      vars &var)
{	
	int n = var.n;
	double err = var.err;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> &urdist = var.urdist;
	boost::random::uniform_real_distribution<> &urdist1 = var.urdist1; 
	
	CGAL::Random_points_on_sphere_d<Point> gen (n, 1.0);
	Vector l = *gen - CGAL::Origin();
	Vector b1 = line_bisect(p,l,P,var);
	Vector b2 = line_bisect(p,-l,P,var);
	//std::cout<<"b1="<<b1<<"b2="<<b2<<std::endl;
	double lambda = urdist(rng);
	p = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
	return 1;
}
*/

// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <class T, class K>
int rand_point_generator(T &P,
												 Point &p,   // a point to start
												 int rnum,
												 int walk_len,
												 K &randPoints,
												 vars &var  // constans for volume
												)
{
  int n = var.n;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> urdist = var.urdist;
	boost::random::uniform_int_distribution<> uidist(0,n-1);
	
	std::vector<NT> lamdas(P.num_of_hyperplanes(),NT(0));
	int rand_coord = uidist(rng);
	double kapa = urdist(rng);
	Point p_prev = p;
	//hit_and_run(p,P,var,var);
	hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord,kapa,lamdas,var,var,true);	
	
  for(int i=1; i<=rnum; ++i){
		
		for(int j=0; j<walk_len; ++j){
		  int rand_coord_prev = rand_coord;
		  rand_coord = uidist(rng);
		  kapa = urdist(rng);
		  //hit_and_run(p,P,var,var);
		  hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord_prev,kapa,lamdas,var,var,false);
		}
		randPoints.push_back(p);		
	}
	//if(rand_only) std::cout<<p<<std::endl;
	//if(print) std::cout<<"("<<i<<") Random point: "<<p<<std::endl;
}

/*
template <class T, class K>
int rand_point_generator_with_walk_estimator(T &P,
														//								 Point &p,   // a point to start
																						 int rnum,
					//																	 int walk_len,
																						 K &randPoints,
																						 vars &var  // constans for volume
																						)
{
  int n = var.n;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> urdist = var.urdist;
	boost::random::uniform_int_distribution<> uidist(0,n-1);
	
	NT walk_estimator=0;
	int walk_len=1;
	while(walk_estimator<0.99){
		Point p=*(randPoints.begin());
		randPoints.clear();
		std::vector<NT> lamdas(P.num_of_hyperplanes(),NT(0));
		int rand_coord = uidist(rng);
		double kapa = urdist(rng);
		Point p_prev = p;
		hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord,kapa,lamdas,var,var,true);	
		
		NT point_coord_sum=0;
		for(int i=1; i<=rnum; ++i){	
			for(int j=0; j<walk_len; ++j){
			  int rand_coord_prev = rand_coord;
			  rand_coord = uidist(rng);
			  kapa = urdist(rng);
			  hit_and_run(p,P,var,var);
			  //hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord_prev,kapa,lamdas,var,var,true);
			}
			randPoints.push_back(p);
			point_coord_sum += p.cartesian(5);
		}
		
		int lpoints=0, rpoints=0; 
		for(typename K::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
			if (pit->cartesian(5) > point_coord_sum/rnum){
				++rpoints;
			} else {
				++lpoints;
			}
		}
		walk_estimator=double(std::min(lpoints,rpoints))/double(std::max(lpoints,rpoints));
		std::cout<<"Walklen="<<walk_len<<" estimator= "<<lpoints<<"/"<<rpoints<<"="<<walk_estimator<<std::endl;
		++walk_len;
	}
	
	//if(rand_only) std::cout<<p<<std::endl;
	//if(print) std::cout<<"("<<i<<") Random point: "<<p<<std::endl;
}
*/

template <class T, class K>
int rand_point_generator(BallIntersectPolytope<T> &PBLarge,
												 Point &p,   // a point to start
												 int rnum,
												 int walk_len,
												 K &randPoints,
												 BallIntersectPolytope<T> &PBSmall,
												 int &nump_PBSmall,
												 vars &var  // constans for volume
												)
{
  int n = var.n;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> urdist = var.urdist;
	boost::random::uniform_int_distribution<> uidist(0,n-1);
	
	std::vector<NT> lamdas(PBLarge.num_of_hyperplanes(),NT(0));
	int rand_coord = uidist(rng);
	double kapa = urdist(rng);
	Point p_prev = p;
	//hit_and_run(p,PBLarge,var,var);
	hit_and_run_coord_update(p,p_prev,PBLarge,rand_coord,rand_coord,kapa,lamdas,var,var,true);	
	
  for(int i=1; i<=rnum; ++i){		
		for(int j=0; j<walk_len; ++j){
		  int rand_coord_prev = rand_coord;
		  rand_coord = uidist(rng);
		  kapa = urdist(rng);
		  //hit_and_run(p,PBLarge,var,var);
		  hit_and_run_coord_update(p,p_prev,PBLarge,rand_coord,rand_coord_prev,kapa,lamdas,var,var,false);
		}
		if(PBSmall.second().is_in(p) == -1){//is in
				randPoints.push_back(p);
				++nump_PBSmall;
		}		
	}
		//if(rand_only) std::cout<<p<<std::endl;
		//if(print) std::cout<<"("<<i<<") Random point: "<<p<<std::endl;
}

// ----- HIT AND RUN FUNCTIONS ------------ //

//hit-and-run with random directions
template <class T>
int hit_and_run(Point &p,
					      T &P,
					      vars &var,
					      vars &var2)
{	
	int n = var.n;
	double err = var.err;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> &urdist = var.urdist;
	boost::random::uniform_real_distribution<> &urdist1 = var.urdist1; 
	
	CGAL::Random_points_on_sphere_d<Point> gen (n, 1.0);
	Vector l = *gen - CGAL::Origin();
	//Vector b1 = line_bisect(p,l,P,var,var2);
	//Vector b2 = line_bisect(p,-l,P,var,var2);
	std::pair<Point,Point> ppair = P.line_intersect(p,l);
	Vector b1 = ppair.first - CGAL::Origin();
	Vector b2 = ppair.second - CGAL::Origin();
	//std::cout<<"b1="<<b1<<"b2="<<b2<<std::endl;
	double lambda = urdist(rng);
	p = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
	return 1;
}

//hit-and-run with orthogonal directions
template <class T>
int hit_and_run_coord(Point &p,
					            T &P,
					            vars &var,
					            vars &var2)
{	
	int n = var.n;
	RNGType &rng = var.rng;
	boost::random::uniform_real_distribution<> &urdist = var.urdist;
	boost::random::uniform_int_distribution<> uidist(0,n-1); 
	
	int rand_coord = uidist(rng);
	double kapa = urdist(rng);
	
	std::pair<NT,NT> bpair = P.line_intersect_coord(p,rand_coord);
	//TODO: only change one coordinate of *r* avoid addition + construction			
	
	std::vector<NT> v(p.dimension(),NT(0));
	v[rand_coord] = bpair.first + kapa * (bpair.second - bpair.first);
	Vector vp(p.dimension(),v.begin(),v.end());
	p = p + vp;

	return 1;
}

//hit-and-run with orthogonal directions and update
template <class T>
int hit_and_run_coord_update(Point &p,
									           Point &p_prev,
					                   T &P,
					                   int rand_coord,
					                   int rand_coord_prev,
					                   double kapa,
					                   std::vector<NT> &lamdas,
					                   vars &var,
					                   vars &var2,
					                   bool init)
{	
	std::pair<NT,NT> bpair;
  if(var.NN)
	  bpair = P.query_dual(p,rand_coord);	
	else
	  bpair = P.line_intersect_coord(p,p_prev,rand_coord,rand_coord_prev,lamdas,init);
	//std::cout<<"original:"<<bpair.first<<" "<<bpair.second<<std::endl;
	//std::cout<<"-----------"<<std::endl;
	//TODO: only change one coordinate of *r* avoid addition + construction			
	std::vector<NT> v(P.dimension(),NT(0));
	v[rand_coord] = bpair.first + kapa * (bpair.second - bpair.first);
	Vector vp(P.dimension(),v.begin(),v.end());
	p_prev = p;
	p = p + vp;
	return 1;
}


/*---------------- MULTIPOINT RANDOM WALK -----------------*/
// generate m random points uniformly distributed in P
template <class T>
int multipoint_random_walk(T &P,
													 std::vector<Point> &V,
											     vars &var)
{
	 int m = var.m;
	 int n = var.n;
	 const int walk_steps = var.walk_steps;
	 const double err = var.err;
	 RNGType &rng = var.rng;
	 generator
	 &get_snd_rand = var.get_snd_rand;
	 boost::random::uniform_real_distribution<> &urdist = var.urdist;
	 boost::random::uniform_real_distribution<> &urdist1 = var.urdist1;
	 
	//remove half of the old points
	//V.erase(V.end()-(V.size()/2),V.end());										
	
	//generate more points (using points in V) in order to have m in total
	std::vector<Point> U;
	std::vector<Point>::iterator Vit=V.begin();
	for(int mk=0; mk<m-V.size(); ++mk){
		// Compute a point as a random uniform convex combination of V 
		//std::vector<double> a;
		//double suma=0;
		//for(int ai=0; ai<V.size(); ++ai){
	  //  a.push_back(urdist(rng));
	  //	suma+=a[a.size()-1];
		//}		
		
		// hit and run at every point in V
		Vector p(n,CGAL::NULL_VECTOR);
	  Point v=*Vit;
	  hit_and_run(v,P,var);
	  U.push_back(v);
	  ++Vit;
	  if(Vit==V.end())
			Vit=V.begin();
	}	
	
	//append U to V
	V.insert(V.end(),U.begin(),U.end());
	//std::cout<<"--------------------------"<<std::endl;
	//std::cout<<"Random points before walk"<<std::endl;
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit){
		Point v=*vit;
		hit_and_run(v,P,var);
		//std::cout<<*vit<<"---->"<<v<<std::endl;
	}
	
	//std::cout<<"WALKING......"<<std::endl;											 
	for(int mk=0; mk<walk_steps; ++mk){
		for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit){
			
			Point v=*vit;
		  
		  // Choose a direction 
			std::vector<double> a(V.size());
			generate(a.begin(),a.end(),get_snd_rand);
			
			std::vector<Point>::iterator Vit=V.begin();
			Vector l(n,CGAL::NULL_VECTOR);
			for(std::vector<double>::iterator ait=a.begin(); ait!=a.end(); ++ait){
			  //*Vit*=*ait;
			  //std::cout<<*ait<<"*"<<(*Vit)<<"= "<<NT(*ait)*(*Vit)<<std::endl;
			  //std::cout<<*ait<<std::endl;
			  l+=NT(*ait)*((*Vit)-(CGAL::Origin()));
			  ++Vit;
			}
			
			// Compute the line 
			Line line(v,l.direction());
			//std::cout<<line<<std::endl;
	    
			// Compute the 2 points that the line and P intersect 
			Vector b1=line_bisect(v,l,P,var);	
			Vector b2=line_bisect(v,-l,P,var);
			//std::cout<<"["<<b1<<","<<b2<<"]"<<std::endl;
			
			// Move the point to a random (uniformly) point in P along the constructed line 
			double lambda = urdist(rng);		
			v = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
			//std::cout<<"new point"<<v<<std::endl;
			//round_print(v);
			*vit=v;
			//hit_and_run(*vit,P,var);
	  }
	}
	/*
	std::cout<<"Random points after walk"<<std::endl;
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit)
		std::cout<<*vit<<std::endl;											 
	std::cout<<"--------------------------"<<std::endl;
	*/
	//for(Polytope::iterator polyit=P.begin(); polyit!=P.end(); ++polyit)
	//	std::cout<<*polyit<<std::endl;
	
	if(m!=V.size()){
		std::cout<<"Careful m!=V.size()!!"<<std::endl;
		exit(1);
	}
}

#endif //RANDOM_SAMPLERS_H
