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


#include <CGAL/point_generators_d.h>
//#include <CGAL/Filtered_kernel_d.h>
//#include <CGAL/Triangulation.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/algorithm.h>
//#include <CGAL/Random.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <bitset>
#include <random>
#include <functional>
#include <algorithm>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"  
#include "boost/dynamic_bitset.hpp"   
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>



//#include <gmpxx.h>
//typedef mpq_class NT;
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
//typedef CGAL::Gmpq                NT;
typedef double                      NT;
//typedef CGAL::Gmpz                NT;


typedef CGAL::Cartesian_d<NT> 	      Kernel; 
//typedef CGAL::Triangulation<Kernel> T;
typedef Kernel::Point_d								Point;
typedef Kernel::Vector_d							Vector;
typedef Kernel::Line_d								Line;
typedef Kernel::Hyperplane_d					Hyperplane;
typedef Kernel::Direction_d						Direction;
//typedef Kernel::Sphere_d						Ball;




// define random generators
typedef boost::mt19937 RNGType; ///< mersenne twister generator
//typedef boost::lagged_fibonacci607 RNGType;

typedef boost::variate_generator< RNGType, boost::normal_distribution<> >  generator;
//typedef boost::variate_generator< RNGType, boost::exponential_distribution<> >  generator;

//structs with variables and random generators
struct vars{
	public:
	  vars( int m,
					int n,
					int walk_steps,
					const double err,
					const double err_opt,
					const int lw,
					double up,
					const int L,
				  RNGType &rng,
				  generator
				  &get_snd_rand,
				  boost::random::uniform_real_distribution<> urdist,
				  boost::random::uniform_real_distribution<> urdist1
			) : 
	    m(m), n(n), walk_steps(walk_steps), err(err), err_opt(err_opt), 
	    lw(lw), up(up), L(L), rng(rng), get_snd_rand(get_snd_rand),
	    urdist(urdist), urdist1(urdist1) {};
	
	  int m;
		int n;
		int walk_steps;
		const double err;
		const double err_opt;
		const int lw;
		double up;
		const int L;
	  RNGType &rng;
	  generator
	  &get_snd_rand;
	  boost::random::uniform_real_distribution<> urdist;
	  boost::random::uniform_real_distribution<> urdist1;
};

// define extreme points
typedef CGAL::Extreme_points_traits_d<Point>   EP_Traits_d;

template <class T>
int optimization(T &KK,vars var,Point &fp,Vector &w);
template <class T>
int opt_interior(T &K,vars &var,Point &opt,Vector &w);

#include <polytopes.h>
#include <oracles.h>
#include <random_samplers.h>
#include <misc.h>





 


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
			Vector b1=line_intersect(v,l,P,var);	
			Vector b2=line_intersect(v,-l,P,var);
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

// return 1 if P is feasible and fp a point in P
// otherwise return 0 and fp has no meaning
template <class T>
int feasibility(T &KK,
								Point &fp,
								vars var)
{
	bool print = false;
	bool print2 = false;
	const int m = var.m;
  int n = var.n;
  const int walk_steps = var.walk_steps;
  const double err = var.err;
  const int lw = var.lw;
  double up = var.up;
  const int L = var.L;
  RNGType &rng = var.rng;
  generator
  &get_snd_rand = var.get_snd_rand;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;
  	
	//this is the large cube contains the polytope
	Polytope P=cube(n,-up,up);								
	
	/* Initialize points in cube */
  // create a vector V with the random points
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(urdist(rng) * up));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		//std::cout<<v<<std::endl;
	}
	
	int step=0;
  while(step < 2*n*L){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,var);
		
	  //compute the average using the half of the random points
		Vector z(n,CGAL::NULL_VECTOR);
		int i=0;
		std::vector<Point>::iterator vit=V.begin();
		//std::cout<<"RANDOM POINTS"<<std::endl;
		for(; i<m/2; ++i,++vit){
			CGAL:assert(vit!=V.end());
			z = z + (*vit - CGAL::Origin());
			//std::cout<<*vit<<std::endl;
		}
		z=z/(m/2);	
		
		if (print) 
		  std::cout<<"Step:"<<step<<" Current centroid= "<<z
		         <<std::endl;
		
		sep sep_result = Sep_Oracle(KK,CGAL::Origin()+z,var);
		if(sep_result.get_is_in()){
			if (print2) std::cout<<"\nFeasible point found! "<<z<<std::endl;
			fp = CGAL::Origin() + z;
			return 1;
		}
		else {
			//update P with the hyperplane passing through z
			Hyperplane H(CGAL::Origin()+z,sep_result.get_H_sep().orthogonal_direction());
			P.push_back(H);
			//GREEDY alternative: Update P with the original separating hyperplane
			//PROBLEM: we may lose all random points thus not efficient
			//Hyperplane H(sep_result.get_H_sep());
			//P.push_back(H);
			
			//check for the rest rand points which fall in new P
			std::vector<Point> newV;
			for(;vit!=V.end();++vit){
				if(Sep_Oracle(P,*vit,var).get_is_in())
					newV.push_back(*vit);
			}
			V=newV;
			++step;
			if (print) std::cout<<"Cutting hyperplane direction="
			         <<sep_result.get_H_sep().orthogonal_direction()<<std::endl;
			if (print) std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
			
			if(V.empty()){//HERE we have not theoretical guarantees
				if (print) std::cout<<"No random points left!"<<std::endl;
				//try to generate new rand points
				for(int i=0; i<m; ++i){
					Point newv = CGAL::Origin() + z;
					hit_and_run(newv,P,var);
					V.push_back(newv);
				}
				//fp = CGAL::Origin() + z;
			}
		}
	}
	if (print2) std::cout<<"\nNo feasible point found!"<<std::endl;
	return 0;
}

/* Optimization with bisection
	 * TODO: make it a function!!!
	 */
 
template <class T>
int opt_bisect(T &K,
							 vars var,
							 Point &opt,
							 Vector &w) 
{
	bool print = false;
	bool print2 = false;
	int m = var.m;
	int n = var.n;
	int walk_steps = var.walk_steps;
	const double err = var.err;
	const double err_opt = var.err_opt;
	const int lw = var.lw;
	double up = var.up;
	const int L = var.L;
  RNGType &rng = var.rng;
  generator
  &get_snd_rand = var.get_snd_rand;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;
	
	
	//compute 2 parallel hyperplanes one feasible and one not
	std::vector<NT> xxin(n,1);
  Point pout(n,xxin.begin(),xxin.end());
  Hyperplane H_out(pout,w);
  //
  std::vector<NT> xxout(n,-1);
  Point pin(n,xxout.begin(),xxout.end());
  Hyperplane H_in(pin,w);
  Point fp(n,xxout.begin(),xxout.end());
  
  //binary search for optimization
  int step = 0;
  double len;
  Point pmid;
  do{
		pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
		Hyperplane H(pmid,w);
		K.push_back(H);
		if (print) std::cout<<"pmid,pin,pout,w"<<std::endl;
		//if (print) round_print(pmid);
		//if (print) round_print(pin);round_print(pout);round_print(w);
		
		NT c0 = w*(fp-CGAL::Origin());
		if(feasibility(K,fp,var) == 1){
			pin=pmid;
		}
		else{
			pout=pmid;
			if (print2) std::cout<<"NON FEASIBLE"<<std::endl;
		}
		K.pop_back();
		
		len=std::abs((pin-CGAL::Origin())*w - (pout-CGAL::Origin())*w);
		//std::cout<<"len="<<len<<std::endl;
		//std::cout<<"fp=";round_print(fp);
	  opt=fp;
    if (print2) std::cout<<"O1 Step:"<<step<<" Current fp= "<<fp
		         <<" w*fp="<<w*(fp-CGAL::Origin())<<" len="<<len<<std::endl;
    step++;
	}while(len > err_opt);
	
	return 1;
}

// TODO: merge same code with interior opt
template <class T>
int optimization(T &KK,
							  vars var,
							  Point &opt,
							  Vector &w)
{
	bool print = false;
	bool print2 = false;
  int m = var.m;
  //std::cout<<"OPT #randpoints="<<m<<std::endl;
	int n = var.n;
	int walk_steps = var.walk_steps;
	const double err = var.err;
	const double err_opt = var.err_opt;
	const int lw = var.lw;
	double up = var.up;
	const int L = var.L;
  RNGType &rng = var.rng;
  generator
  &get_snd_rand = var.get_snd_rand;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

	//this is the bounding cube that contains the polytope KK
	Polytope P=cube(n,-up,up);							
	
	/* Initialize points in cube */  
  // create a vector V with the random points
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(urdist(rng) * up));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		//std::cout<<v<<std::endl;
	}
	
	//initialize the cut with sth that contain KK
	Hyperplane KK_cut = *(P.begin());
	//iterate for 2nL steps 
  int step=0;
  double epsilon=1.0;
  Vector z(n,CGAL::NULL_VECTOR);
    
  while(step < 2*n*L && epsilon > err_opt){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,var);
		
	
	  //compute the average using the half of the random points
		Vector newz(n,CGAL::NULL_VECTOR);
		int i=0;
		std::vector<Point>::iterator vit=V.begin();
		if (print) std::cout<<"RANDOM POINTS"<<std::endl;
		for(; i<m/2; ++i,++vit){
			CGAL:assert(vit!=V.end());
			newz = newz + (*vit - CGAL::Origin());
			//std::cout<<*vit<<std::endl;
		}
		
		newz=newz/(m/2);	
		epsilon = std::abs(w*newz - w*z);
		
		//update z
		z = newz;
		
		//std::cout<<"step "<<step<<": "<<"z=";
		if (print) round_print(z);
		if (print2) std::cout<<"O2 Step:"<<step<<" Current centroid= "<<z
		         <<" w*z="<<w*z<<" epsilon="<<epsilon<<std::endl;
		
		sep sep_result = Sep_Oracle(KK,CGAL::Origin()+z,var);
		if(sep_result.get_is_in()){
			if (print) std::cout<<"Feasible point found! "<<z<<std::endl;
			opt = CGAL::Origin() + z;
			Hyperplane H(opt,w);
			P.push_back(H);
			//KK.push_back(H);
		}
		else {
			//update P with the hyperplane passing through z
			Hyperplane H(CGAL::Origin()+z,sep_result.get_H_sep().orthogonal_direction());
			P.push_back(H);
			//GREEDY alternative: Update P with the original separating hyperplane
			//Hyperplane H(sep_result.get_H_sep());
			//P.push_back(H);
		}	
		//check for the rest rand points which fall in new P
		std::vector<Point> newV;
		for(;vit!=V.end();++vit){
			if(Sep_Oracle(P,*vit,var).get_is_in())
				newV.push_back(*vit);
		}
		V=newV;
		++step;
		if (print) std::cout<<"Cutting hyperplane direction="
		         <<sep_result.get_H_sep().orthogonal_direction()<<std::endl;
		if (print) std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
		
		if(V.empty()){//HERE we have no theoretical guarantees
				if (print) std::cout<<"No random points left!"<<std::endl;
				//try to generate new rand points
				for(int i=0; i<m; ++i){
					Point newv = CGAL::Origin() + z;
					hit_and_run(newv,P,var);
					V.push_back(newv);
				}
			//opt = CGAL::Origin() + z;
		}
	}
	return 0;
}

// TODO: merge same code with interior opt
template <class T>
int optimization2(T &KK,
							    vars var,
							    Point &opt,
							    Vector &w)
{
	bool print = false;
	bool print2 = false;
  int m = var.m;
  //std::cout<<"OPT #randpoints="<<m<<std::endl;
	int n = var.n;
	int walk_steps = var.walk_steps;
	const double err = var.err;
	const double err_opt = var.err_opt;
	const int lw = var.lw;
	double up = var.up;
	const int L = var.L;
  RNGType &rng = var.rng;
  generator
  &get_snd_rand = var.get_snd_rand;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

	//this is the bounding cube that contains the polytope KK
	Polytope P=cube(n,-up,up);							
	
	/* Initialize points in cube */  
  // create a vector V with the random points
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(urdist(rng) * up));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		//std::cout<<v<<std::endl;
	}
	
	//compute 2 parallel hyperplanes one feasible and one not
	std::vector<NT> xxin(n,up);
  Point pout(n,xxin.begin(),xxin.end());
  Hyperplane H_out(pout,w);
  //
  std::vector<NT> xxout(n,-up);
  Point pin(n,xxout.begin(),xxout.end());
  Hyperplane H_in(pin,w);
  //
  Point pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
  Hyperplane H_mid(pmid,w);
  
	//initialize the cut with sth that contain KK
	Hyperplane KK_cut = *(P.begin());
	//iterate for 2nL steps 
  int step=0, small_step=0;
  double epsilon=1.0;
  Vector z(n,CGAL::NULL_VECTOR);
    
  //while(step < 2*n*L && epsilon > err_opt){
	while(epsilon > err_opt){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,var);
		
	  //compute the average using the half of the random points
		Vector newz(n,CGAL::NULL_VECTOR);
		int i=0;
		std::vector<Point>::iterator vit=V.begin();
		if (print) std::cout<<"RANDOM POINTS"<<std::endl;
		for(; i<m/2; ++i,++vit){
			CGAL:assert(vit!=V.end());
			newz = newz + (*vit - CGAL::Origin());
			//std::cout<<*vit<<std::endl;
		}
		
		newz=newz/(m/2);	
		epsilon = std::abs(w*(pin-CGAL::Origin()) - w*(pout-CGAL::Origin()));
		
		//update z
		z = newz;
		
		//std::cout<<"step "<<step<<": "<<"z=";
		if (print) round_print(z);
		if (print2) std::cout<<"Step:"<<step<<" Current centroid= "<<z
		         <<" w*z="<<w*z<<" epsilon="<<epsilon<<std::endl;
		
		sep sep_result = Sep_Oracle(KK,CGAL::Origin()+z,var);
		if(sep_result.get_is_in()){
			if (print) std::cout<<"Feasible point found! "<<z<<std::endl;
			opt = CGAL::Origin() + z;
			pin=opt;
			H_in = Hyperplane(pin,w);
			P.push_back(H_in);
			//KK.push_back(H);
			if (!H_mid.has_on_negative_side(CGAL::Origin()+z)){
				pin=pmid;
				H_in = Hyperplane(pin,w);
				pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
			  H_mid = Hyperplane(pmid,w);
			}
			pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
			H_mid = Hyperplane(pmid,w);
		}
		else{
			//update P with the hyperplane passing through z
			Hyperplane H(CGAL::Origin()+z,sep_result.get_H_sep().orthogonal_direction());
			P.push_back(H);
			//GREEDY alternative: Update P with the original separating hyperplane
			//Hyperplane H(sep_result.get_H_sep());
			//P.push_back(H);
		}
		if(small_step > 2*n){
			//H_mid is infeasible. update it
			pout=pmid;
			H_out = Hyperplane(pout,w);
			pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
			H_mid = Hyperplane(pmid,w);
			small_step=0;
		}	
		//check for the rest rand points which fall in new P
		std::vector<Point> newV;
		for(;vit!=V.end();++vit){
			if(Sep_Oracle(P,*vit,var).get_is_in())
				newV.push_back(*vit);
		}
		V=newV;
		++step;
		++small_step;
		if (print) std::cout<<"Cutting hyperplane direction="
		         <<sep_result.get_H_sep().orthogonal_direction()<<std::endl;
		if (print) std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
		
		if(V.empty()){//HERE we have no theoretical guarantees
				if (print) std::cout<<"No random points left!"<<std::endl;
				//try to generate new rand points
				for(int i=0; i<m; ++i){
					Point newv = CGAL::Origin() + z;
					hit_and_run(newv,P,var);
					V.push_back(newv);
				}
			//opt = CGAL::Origin() + z;
		}
	}
	return 0;
}

// interior point optimization
template <class T>
int opt_interior(T &K,
                 vars &var,
							   Point &opt,
							   Vector &w){
	
	bool print = false;
	bool print2 = false;
	int m = var.m;
	int n = var.n;
	int walk_steps = var.walk_steps;
	const double err = var.err;
	const double err_opt = var.err_opt;
	const int lw = var.lw;
	double up = var.up;
	const int L = var.L;
  RNGType &rng = var.rng;
  generator
  &get_snd_rand = var.get_snd_rand;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;
	
	//first compute a feasible point in K
	Point fp;
  if (feasibility(K,fp,var)==0){
	  if (print) std::cout<<"The input polytope is not feasible!"<<std::endl;
	  return 1;
	}
	Vector z = fp - CGAL::Origin();
	// Initialization
	Hyperplane H(fp,w);
	K.push_back(H);
	
	// create a vector V with the random points
	std::vector<Point> V;
	//initialize V !!!! This is a case without theoretical guarantees 
	for(int i=0; i<m; ++i){
		Point newv = CGAL::Origin() + z;
		hit_and_run(newv,K,var);
	  V.push_back(newv);
	}
	//
	int step=0;
	double epsilon;	
	do{
		
		// compute m random points in K stored in V 
	  multipoint_random_walk(K,V,var);
			
	  //compute the average (new z) using the half of the random points
		Vector newz(n,CGAL::NULL_VECTOR);
		
		int i=0;
		std::vector<Point>::iterator vit=V.begin();
		//std::cout<<"RANDOM POINTS"<<std::endl;
		for(; i<m/2; ++i,++vit){
			CGAL:assert(vit!=V.end());
			newz = newz + (*vit - CGAL::Origin());
			//std::cout<<*vit<<std::endl;
		}
		
		newz=newz/(m/2);	
		epsilon = std::abs(w*newz - w*z);
		
		if (print2) 
		  std::cout<<"O3 Step:"<<step<<" Current centroid= "<<z
		         <<" w*z="<<w*z<<" epsilon="<<epsilon<<std::endl;
		
		//Update z
		z=newz;
		
		//update K with the hyperplane passing through z	
		Hyperplane H(CGAL::Origin()+z,w);
		K.pop_back();
		K.push_back(H);
		
		//check for the rest rand points which fall in new K
		std::vector<Point> newV;
		for(;vit!=V.end();++vit){
			if(Sep_Oracle(K,*vit,var).get_is_in())
				newV.push_back(*vit);
		}
		V=newV;
		++step;
		
		if (print) { 
	    std::cout<<"Cutting hyperplane direction="
			         <<H.orthogonal_direction()<<std::endl;
		  std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
	  }
		if(V.empty()){
			//HERE we have no theoretical guarantees
			if (print) std::cout<<"No random points left!"<<std::endl;
			//try to generate new rand points
			for(int i=0; i<m; ++i){
				Point newv = CGAL::Origin() + z;
				hit_and_run(newv,K,var);
				V.push_back(newv);
			}
		}
		
	}while(epsilon>err_opt);
	
	if (print) std::cout<<"OPT = "<<z<<std::endl;
	opt = CGAL::Origin() + z;
	return 0;
}


/////////////////////////////////////////////////////////
// VOLUME
// randomized approximate volume computation 
/*************************************************
/**************** VOLUME with hit and run        */
// We assume that the polytope P is properly sandwitched
// The sandwitching:
// r is the radius of the smallest ball
// d is the radius of the largest
template <class T>
NT volume1(T &P,
					 vars &var,  // constans for volume
					 vars &var2, // constants for optimization in case of MinkSums
					 NT r, 
					 NT d)
{
  typedef BallIntersectPolytope<T>        BallPoly; 
	
	bool print = true;
	int n = var.n;
	int rnum = var.m;
	int walk_len = var.walk_steps;
	const double err = var.err;
	RNGType &rng = var.rng;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;
	
	//The number of balls we construct 
  int nb = std::ceil(n * (std::log(d)/std::log(2.0)));
  //std::cout<<"nb="<<nb<<", d="<<d<<std::endl;
  //std::pow(std::log(n),2)
  
  // Construct the sequence of balls
  std::vector<NT> coords(n,0);
	Point p0(n,coords.begin(),coords.end());
  std::vector<Ball> balls;
  for(int i=0; i<=nb; ++i){
		balls.push_back(Ball(p0,std::pow(std::pow(2.0,NT(i)/NT(n)),2))); 
		if (print) std::cout<<"ball"<<i<<"="<<balls[i].center()<<" "<<balls[i].radius()<<std::endl;
	}
  assert(!balls.empty());
  if (print) std::cout<<"---------"<<std::endl;
  
    
  std::vector<int> telescopic_prod(nb,0);
  for(int i=1; i<=rnum; ++i){ //generate rnum rand points 
		//start with a u.d.r point in the smallest ball 
		//radius=1, center=Origin()
		std::vector<NT> coords(n,0);
		Point p(n,coords.begin(),coords.end());
		BallPoly PBold(P,balls[0]);
		//
		//hit_and_run(p,PBold.second(),var,var2);
		CGAL::Random_points_in_ball_d<Point> gen (n, 1.0);
		p = *gen;
		//std::cout<<p<<std::endl;
		//std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
		//std::cout<<balls[0].is_in(p)<<std::endl;
		
		//exit(0);
		std::vector<Ball>::iterator bit=balls.begin();
		std::vector<int>::iterator prod_it=telescopic_prod.begin();
		++bit; 
		for(; bit!=balls.end(); ++bit, ++prod_it){
			// generate a random point in bit intersection with P 
			BallPoly PB(P,*bit);
			
			for(int j=0; j<walk_len; ++j){
			  hit_and_run(p,PB,var,var2);
				//std::cout<<"h-n-r:"<<p<<std::endl;
			}
			//Not need to test for PBold membership. Just check if inside Ball
			//if (Sep_Oracle(PBold,p,var2).get_is_in()){
			if (PBold.second().is_in(p)){
				//std::cout<<p<<" IN ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
			  ++(*prod_it);
			}else{
			  ;
			  //std::cout<<p<<":"<<(p-CGAL::Origin()).squared_length()
			  //<<" OUT ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
			}
			PBold=PB;
		}
		if (print) std::cout<<"\n\ngenerated random point..."<<i<<"/"<<rnum<<" ";
		const NT pi = boost::math::constants::pi<NT>();
		NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0); 
		for(std::vector<int>::iterator prod_it=telescopic_prod.begin(); 
		    prod_it!=telescopic_prod.end(); ++prod_it){
			vol *= NT(i)/NT(*prod_it);
		}
		if (print) std::cout<<"current vol estimation= "<<vol<<std::endl;
	  if (print) std::cout<<"walklen="<<walk_len<<std::endl;
	  if (print) std::cout<<"rnum="<<rnum<<std::endl;
		//for(prod_it=telescopic_prod.begin(); prod_it!=telescopic_prod.end(); ++prod_it)
		//	std::cout<<(*prod_it)<<" ";
		//std::cout<<std::endl;
	}	
	const NT pi = boost::math::constants::pi<NT>();
	NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0); 
	//NT vol=1;
	if (print) std::cout<<"vol(K_0)="<<vol<<" ";
	for(std::vector<int>::iterator prod_it=telescopic_prod.begin();
	    prod_it!=telescopic_prod.end(); ++prod_it){
		vol *= NT(rnum)/NT(*prod_it);
		if (print) std::cout<<NT(rnum)<<"/" << NT(*prod_it)<<"="<<NT(rnum)/NT(*prod_it)<<"\n";
	}
	return vol;
}


// VOLUME with multipoint random walk
template <class T>
NT volume2(T &P,
					 vars &var)
{
	typedef BallIntersectPolytope<T>        BallPoly; 				
				
	int n = var.n;
	int rnum = var.m;
	int walk_len = var.walk_steps;
	const double err = var.err;
	RNGType &rng = var.rng;
  generator
  get_snd_rand = var.get_snd_rand;
  boost::random::uniform_real_distribution<> urdist = var.urdist;
  boost::random::uniform_real_distribution<> urdist1 = var.urdist1;
  	
	// The sandwitching 
	// r is the radius of the smallest ball
	// d is the radius of the largest
  std::vector<NT> coords_apex(n,1);
	Vector p_apex(n,coords_apex.begin(),coords_apex.end());
  const NT r=1, d=std::sqrt(p_apex.squared_length());
  const int nb = std::ceil(n * (std::log(d)/std::log(2.0)));
  std::cout<<"nb="<<nb<<", d="<<d<<std::endl;
  //std::pow(std::log(n),2)
  
  // Construct the sequence of balls
  std::vector<NT> coords(n,0);
	Point p0(n,coords.begin(),coords.end());
  std::vector<Ball> balls;
  for(int i=0; i<=nb; ++i){
		balls.push_back(Ball(p0,std::pow(std::pow(2.0,NT(i)/NT(n)),2))); 
		std::cout<<"ball"<<i<<"="<<balls[i].center()<<" "<<balls[i].radius()<<std::endl;
	}
  assert(!balls.empty());
  std::cout<<"---------"<<std::endl;
  
  
  
  std::vector<int> telescopic_prod(nb,0);
  //vector to store the random points
  std::vector<Point> V;
  BallPoly PBold(P,balls[0]);
  for(int i=0; i<rnum; ++i){ 
		// generate rnum rand points  
		// in the smallest ball i.e. radius=1, center=Origin()
		std::vector<NT> coords(n,0);
		Point p(n,coords.begin(),coords.end());
		hit_and_run(p,PBold,var);
		V.push_back(p);
		//std::cout<<p<<std::endl;
	}	
	//std::cout<<p<<std::endl;
	//std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
	//std::cout<<balls[0].is_in(p)<<std::endl;
	//exit(0);
	std::vector<Ball>::iterator bit=balls.begin();
	std::vector<int>::iterator prod_it=telescopic_prod.begin();
	++bit; 
	for(; bit!=balls.end(); ++bit, ++prod_it){
		// generate a random point in bit (intersection) P 
		BallPoly PB(P,*bit);
	  std::cout<<"walking..."<<walk_len<<"steps"<<std::endl;
	  var.m = V.size();
		multipoint_random_walk(PB,V,var);
		
		for(int j=0; j<V.size(); ++j){
		  
			if (Sep_Oracle(PBold,V[j],var).get_is_in()){
				//std::cout<<p<<" IN ball: "<<PBold.second.center()<<PBold.second.radius()<<std::endl;
				++(*prod_it);
			}//else{
				//std::cout<<p<<":"<<(p-CGAL::Origin()).squared_length()
				//<<" OUT ball: "<<PBold.second.center()<<PBold.second.radius()<<std::endl;
		}
		PBold=PB;	
		//for(prod_it=telescopic_prod.begin(); prod_it!=telescopic_prod.end(); ++prod_it)
		//	std::cout<<(*prod_it)<<" ";
		//std::cout<<std::endl;
	}	
	const NT pi = boost::math::constants::pi<NT>();
	NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0); 
	//NT vol=1;
	std::cout<<"vol(K_0)="<<vol<<" ";
	for(std::vector<int>::iterator prod_it=telescopic_prod.begin(); 
	    prod_it!=telescopic_prod.end(); ++prod_it){
		vol *= NT(rnum)/NT(*prod_it);
		std::cout<<NT(rnum)<<"/" << NT(*prod_it)<<"="<<NT(rnum)/NT(*prod_it)<<"\n";
	}
	std::cout<<std::endl;
	std::cout<<"volume = "<<vol<<std::endl;
	std::cout<<"exact volume = "<<std::pow(2,n)<<std::endl;
	return vol;
}


// polymake file to compute exact volume
template <class T>
void print_polymake_volfile(T &P,
														std::ostream& os){
  // print the vertices of the P polytope
  os << "use Time::HiRes qw(gettimeofday tv_interval);\n";
	os << "use application 'polytope';\n";
	os << "my $p=new Polytope<Rational>;\n";
	os << "$p->POINTS=<<'.';\n";
  for (typename T::iterator vit = P.begin(); vit != P.end(); vit++){
    os << "1 ";
    for (Point::Cartesian_const_iterator cit=vit->cartesian_begin();
         cit != vit->cartesian_end();
         cit++){
      os << *cit;
      if (cit - vit->cartesian_begin() != vit->dimension()-1)
        os << " ";
    }
    //os << "|" << vit->point().index();
    os << "\n";
  }
	os << ".\n";
	os << "print ' ';\n";
	os << "print $p->N_POINTS;\n";
	os << "print ' ';\n";
	os << "print $p->N_VERTICES;\n";
	os << "print ' ';\n";
	os << "print $p->DIM;\n";
	os << "print ' ';\n";
	os << "my $t0 = [gettimeofday];\n";
	os << "my $f=$p->VOLUME;\n";
	os << "print $f;\n";
	os << "print ' ';\n";
	os << "print tv_interval($t0,[gettimeofday]);\n";
	os << "print \"\n\";\n";
	os << std::endl;
}

// polymake file to compute exact volume
template <class T>
void print_polymake_volfile2(T &P,
														std::ostream& os){
  // print the vertices of the P polytope
  os << "use Time::HiRes qw(gettimeofday tv_interval);\n";
	os << "use application 'polytope';\n";
	os << "my $p=new Polytope;\n";
	os << "$p->INEQUALITIES=<<'.';\n";
	//os << "my $p=new Polytope<Rational>;\n";
	//os << "$p->POINTS=<<'.';\n";
  for (typename T::iterator vit = P.begin(); vit != P.end(); vit++){
    Hyperplane::Coefficient_const_iterator cit_end = vit->coefficients_end();
    os << *(--cit_end)<<" ";
    //os << "0 ";
    Hyperplane::Coefficient_const_iterator cit = vit->coefficients_begin();
    //++cit;
    for (; cit != cit_end; cit++){
      //std::cout<<*cit<<" ";
      os <<(*cit)<<" ";
      if (cit - vit->coefficients_begin() != vit->dimension()-1)
        os << " ";
    }
    //os << "|" << vit->point().index();
    os << "\n";
  }
	os << ".\n";
	//$p=new Polytope<Rational>(INEQUALITIES=>$inequalities);

	os << "print ' ';\n";
	os << "print $p->N_POINTS;\n";
	os << "print ' ';\n";
	os << "print $p->N_VERTICES;\n";
	os << "print ' ';\n";
	os << "print $p->DIM;\n";
	os << "print ' ';\n";
	os << "my $t0 = [gettimeofday];\n";
	os << "my $f=$p->VOLUME;\n";
	os << "print $f;\n";
	os << "print ' ';\n";
	os << "print tv_interval($t0,[gettimeofday]);\n";
	os << "print \"\n\";\n";
	os << std::endl;
}

