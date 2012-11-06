#include <CGAL/point_generators_d.h>
//#include <CGAL/Filtered_kernel_d.h>
//#include <CGAL/Triangulation.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/algorithm.h>
//#include <CGAL/Random.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"    
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>

//#include <gmpxx.h>
//typedef mpq_class NT;
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
//typedef CGAL::Gmpq                NT;
typedef double                NT;
//typedef CGAL::Gmpz                NT;


typedef CGAL::Cartesian_d<NT> 	      Kernel; 
//typedef CGAL::Triangulation<Kernel> T;
typedef Kernel::Point_d								Point;
typedef Kernel::Vector_d							Vector;
typedef Kernel::Line_d								Line;
typedef Kernel::Hyperplane_d					Hyperplane;
typedef Kernel::Direction_d						Direction;
typedef std::vector<Hyperplane>       H_polytope;
typedef H_polytope 								    Polytope;
typedef std::vector<Point>						V_polytope;
typedef std::pair<V_polytope,V_polytope> 	MinkSumPolytope;

typedef CGAL::Extreme_points_traits_d<Point>   EP_Traits_d;

typedef boost::mt19937 RNGType; ///< mersenne twister generator

//function to print rounding to double coordinates 
template <class T>
void round_print(T p) { 
  for(typename T::Cartesian_const_iterator cit=p.cartesian_begin(); 
      cit!=p.cartesian_end(); ++cit)
	  std::cout<<CGAL::to_double(*cit)<<" "; 
  std::cout<<std::endl;
}

// Naive algorithm for Mink sum
typedef std::vector<V_polytope>  	Vpolys;

int Minkowski_sum_naive(V_polytope &P1, V_polytope &P2, V_polytope &Msum){
	std::cout<<(!P1.empty() && !P2.empty())<<std::endl;
	if(!P1.empty() && !P2.empty()){
	  V_polytope Msum_all;
		for (V_polytope::iterator Pit1 = P1.begin(); Pit1 != P1.end(); ++Pit1){
	    for (V_polytope::iterator Pit2 = P2.begin(); Pit2 != P2.end(); ++Pit2){
	      Point p = CGAL::Origin() + 
	            (((*Pit1)-CGAL::Origin()) + ((*Pit2)-CGAL::Origin()));
	      Msum_all.push_back(p);
	      std::cout<<p<<std::endl;
	    }
	  } 
	  std::cout<<"---------"<<std::endl;
	  // compute the extreme points
		CGAL::Extreme_points_d<EP_Traits_d> ep(P1[0].dimension());
	  ep.insert(Msum_all.begin(),Msum_all.end());
		//std::vector<Point> extreme_points;
		ep.get_extreme_points(std::back_inserter(Msum));
		return Msum.size();
  }
  return -1;
}

// separation oracle return type 
struct sep{
	public:
	  sep(bool in, Hyperplane H) : is_in(in), H_sep(H) {}
	  sep(bool in) : is_in(in), H_sep(Hyperplane()) {}
	  
	  bool get_is_in(){
			return is_in;
		}
		Hyperplane get_H_sep(){
			return H_sep;
		}
	private:
		bool 				is_in;
		Hyperplane 	H_sep;
};


/* Construct a n-CUBE */
Polytope cube(const int n, const int lw, const int up){	
	Polytope cube;
	std::vector<NT> origin(n,NT(lw));
	for(int i=0; i<n; ++i){
		std::vector<NT> normal;
		for(int j=0; j<n; ++j){
			if(i==j) 
				normal.push_back(NT(1));
			else normal.push_back(NT(0));
		}
		Hyperplane h(Point(n,origin.begin(),origin.end()),
	           Direction(n,normal.begin(),normal.end()));
	  cube.push_back(h);
	}
	std::vector<NT> apex(n,NT(up));
	for(int i=0; i<n; ++i){
		std::vector<NT> normal;
		for(int j=0; j<n; ++j){
			if(i==j) 
				normal.push_back(NT(-1));
			else normal.push_back(NT(0));
		}
		Hyperplane h(Point(n,apex.begin(),apex.end()),
	           Direction(n,normal.begin(),normal.end()));
	  cube.push_back(h);
	}
	return cube;
}

//template <typename T> struct Oracle{
//  sep Sep_Oracle(T &P, Point v);
//};

template <typename T> sep Sep_Oracle(T&, Point);

//function that implements the separation oracle 
template<> sep Sep_Oracle<Polytope>(Polytope &P, Point v)
{
	typename Polytope::iterator Hit=P.begin(); 
	while (Hit!=P.end()){
		if (Hit->has_on_negative_side(v))
			return sep(false,*Hit);
		++Hit;
	}
	return sep(true);	
}
 
// Minkowski sum Separation 
template<> sep Sep_Oracle<MinkSumPolytope>(MinkSumPolytope &P, 
																				Point v)
{		
	Vector q=v-CGAL::Origin();
	V_polytope P1=P.first;
	V_polytope P2=P.second;
	NT max = q * (*(P1.begin())-CGAL::Origin());
	Point max_p1 = *(P1.begin());
	for(V_polytope::iterator pit=P1.begin(); pit!=P1.end(); ++pit){
		double innerp = q * (*pit-CGAL::Origin());
		std::cout<<*pit<<" "<<q<<" "<<innerp<<std::endl;
		if(innerp > max){
			max = innerp;
			max_p1 = *pit;
		}
	}
	std::cout<<max_p1<<std::endl;
	max = q * (*(P2.begin())-CGAL::Origin());
	Point max_p2 = *(P2.begin());
	for(V_polytope::iterator pit=P2.begin(); pit!=P2.end(); ++pit){
		double innerp = q * (*pit-CGAL::Origin());
		if(innerp > max){
			max = innerp;
			max_p2 = *pit;
		}
	}
	std::cout<<max_p2<<std::endl;
	Vector max_psum=(max_p1-CGAL::Origin())+(max_p2-CGAL::Origin());
	//+ q*(-P1sum -P2sum)
	if(max_psum*q <= 1)
		return sep(true);	
	else{
		//construct the dual hyperplane of max_psum
	  std::vector<NT> vcoord;
	  for(Vector::Cartesian_const_iterator cit=max_psum.cartesian_begin(); 
        cit!=max_psum.cartesian_end(); ++cit){
		  vcoord.push_back(*cit);
		}
	  Hyperplane H_sep(int(P1.size()),vcoord.begin(),vcoord.end());
	  return sep(false,H_sep);
	}
	//if(max_psum*q == 1)
	//	std::cout<<"sharp!"<<std::endl;
}	
 
 
// function to find intersection of a line and a polytope 
Vector line_intersect(Point pin, Vector l, Polytope P, double err){
  Vector vin=pin-CGAL::Origin();
  //first compute a point outside P along the line
  Point pout=pin;
  //std::cout<<"Starting inside point: ";
  //round_print(pin);
  
  Vector aug(l);
  while(Sep_Oracle(P,pout).get_is_in() == true){
    aug*=2;
    pout+=aug;
    //std::cout<<"Outside point: ";
    //round_print(pout);
  }
  Vector vout=pout-CGAL::Origin();
  
  //intersect using bisection
  //std::cout<<"pout"<<vout<<std::endl;
  Vector vmid;
  double len;
  do{
		vmid=(vin+vout)/2;
		if(Sep_Oracle(P,CGAL::Origin()+vmid).get_is_in() == false)
			vout=vmid;
		else
			vin=vmid;
		len=CGAL::to_double((vin-vout).squared_length());
		//std::cout<<"len="<<bool(len<err)<<std::endl;
	}while(len > err);
  
  //std::cout<<"Intersection point: ";
  //round_print(vmid);
  return vmid; 
}

/* Hit and run Random Walk */
int hit_and_run(Point &p,
					      Polytope &P,
					      const int n,
					      double err,
								RNGType &rng,
								boost::random::uniform_real_distribution<> &urdist,
								boost::random::uniform_real_distribution<> &urdist1){
									
	
	std::vector<NT> v;
	for(int i=0; i<n; ++i)
		v.push_back(urdist1(rng));
	Vector l(n,v.begin(),v.end());
	Vector b1 = line_intersect(p,l,P,err);
	Vector b2 = line_intersect(p,-l,P,err);
	double lambda = urdist(rng);
	p = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
	return 1;
}

/*---------------- MULTIPOINT RANDOM WALK -----------------*/
// generate m random points uniformly distributed in P
int multipoint_random_walk(Polytope &P,
													 std::vector<Point> &V,
													 const int m,
													 const int n,
													 const int walk_steps,
													 const double err,
													 RNGType &rng,
													 boost::variate_generator< RNGType, boost::normal_distribution<> >
													 &get_snd_rand,
													 boost::random::uniform_real_distribution<> &urdist,
													 boost::random::uniform_real_distribution<> &urdist1){
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
	  hit_and_run(v,P,n,err,rng,urdist,urdist1);
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
		hit_and_run(v,P,n,err,rng,urdist,urdist1);
		//std::cout<<*vit<<"---->"<<v<<std::endl;
	}
	//std::cout<<"WALKING......"<<std::endl;											 
	for(int mk=0; mk<walk_steps; ++mk){
		for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit){
			Point v=*vit;
			
			/* Choose a direction */
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
			Vector b1=line_intersect(v,l,P,err);
			Vector b2=line_intersect(v,-l,P,err);
			//std::cout<<"["<<b1<<","<<b2<<"]"<<std::endl;
			
			// Move the point to a random (uniformly) point in P along the constructed line 
			double lambda = urdist(rng);		
			v = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
			//std::cout<<"new point"<<v<<std::endl;
			//round_print(v);
			*vit=v;
	  }
	}
	/*
	std::cout<<"Random points after walk"<<std::endl;
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit)
		std::cout<<*vit<<std::endl;											 
	std::cout<<"--------------------------"<<std::endl;
	*/
	for(Polytope::iterator polyit=P.begin(); polyit!=P.end(); ++polyit)
		std::cout<<*polyit<<std::endl;
	
	if(m!=V.size()){
		std::cout<<"Careful m!=V.size()!!"<<std::endl;
		exit(1);
	}
}

// return 1 if P is feasible and fp a point in P
// otherwise return 0 and fp has no meaning
int feasibility(Polytope &KK,
							  const int m,
							  const int n,
							  const int walk_steps,
							  const double err,
							  const int lw,
							  const int up,
							  const int L,
							  RNGType &rng,
							  boost::variate_generator< RNGType, boost::normal_distribution<> >
							  &get_snd_rand,
							  boost::random::uniform_real_distribution<> urdist,
							  boost::random::uniform_real_distribution<> urdist1,
							  Point &fp){
	
	//this is the large cube contains the polytope
	Polytope P=cube(n,lw,up);								
	
	/* Initialize points in cube */
  CGAL::Random CGALrng;
  // create a vector V with the random points
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(CGALrng.get_int(lw,up)));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		//std::cout<<v<<std::endl;
	}
	
	int step=0;
  while(step < 2*n*L){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,m,n,walk_steps,err,rng,get_snd_rand,urdist,urdist1);
		
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
		std::cout<<"step "<<step<<": "<<"z=";
		round_print(z);
		
		sep sep_result = Sep_Oracle(KK,CGAL::Origin()+z);
		if(sep_result.get_is_in()){
			std::cout<<"Feasible point found! "<<z<<std::endl;
			fp = CGAL::Origin() + z;
			return 1;
		}
		else {
			//update P with the hyperplane passing through z
			Hyperplane H(CGAL::Origin()+z,sep_result.get_H_sep().orthogonal_direction());
			P.push_back(H);
			//GREEDY alternative: Update P with the original separating hyperplane
			//Hyperplane H(sep_result.get_H_sep());
			//P.push_back(H);
			
			//check for the rest rand points which fall in new P
			std::vector<Point> newV;
			for(;vit!=V.end();++vit){
				if(Sep_Oracle(P,*vit).get_is_in())
					newV.push_back(*vit);
			}
			V=newV;
			++step;
			std::cout<<"Cutting hyperplane direction="
			         <<sep_result.get_H_sep().orthogonal_direction()<<std::endl;
			std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
			if(V.empty()){
				std::cout<<"No random points left. ASSUME that there is no feasible point!"<<std::endl;
				//fp = CGAL::Origin() + z;
				return 0;
			}
		}
	}
	std::cout<<"No feasible point found!"<<std::endl;
	return 0;
}

// return 1 if P is feasible and fp a point in P
// otherwise return 0 and fp has no meaning
int optimization(Polytope &KK,
							  const int m,
							  const int n,
							  const int walk_steps,
							  const double err,
							  const int lw,
							  const int up,
							  const int L,
							  RNGType &rng,
							  boost::variate_generator< RNGType, boost::normal_distribution<> >
							  &get_snd_rand,
							  boost::random::uniform_real_distribution<> urdist,
							  boost::random::uniform_real_distribution<> urdist1,
							  Point &fp,
							  Vector &w){
	
	//this is the large cube contains the polytope
	Polytope P=cube(n,lw,up);								
	
	/* Initialize points in cube */
  CGAL::Random CGALrng;
  // create a vector V with the random points
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(CGALrng.get_int(lw,up)));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		//std::cout<<v<<std::endl;
	}
	
	int step=0;
  while(step < 2*n*L){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,m,n,walk_steps,err,rng,get_snd_rand,urdist,urdist1);
		
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
		std::cout<<"step "<<step<<": "<<"z=";
		round_print(z);
		
		sep sep_result = Sep_Oracle(KK,CGAL::Origin()+z);
		if(sep_result.get_is_in()){
			std::cout<<"Feasible point found! "<<z<<std::endl;
			fp = CGAL::Origin() + z;
			Hyperplane H(fp,w);
			KK.push_back(H);
		}
		else {
			//update P with the hyperplane passing through z
			Hyperplane H(CGAL::Origin()+z,sep_result.get_H_sep().orthogonal_direction());
			P.push_back(H);
			//GREEDY alternative: Update P with the original separating hyperplane
			//Hyperplane H(sep_result.get_H_sep());
			//P.push_back(H);
			
			//check for the rest rand points which fall in new P
			std::vector<Point> newV;
			for(;vit!=V.end();++vit){
				if(Sep_Oracle(P,*vit).get_is_in())
					newV.push_back(*vit);
			}
			V=newV;
			++step;
			std::cout<<"Cutting hyperplane direction="
			         <<sep_result.get_H_sep().orthogonal_direction()<<std::endl;
			std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
			if(V.empty()){
				std::cout<<"No random points left. ASSUME that there is no feasible point!"<<std::endl;
				//fp = CGAL::Origin() + z;
				return 0;
			}
		}
	}
	std::cout<<"No feasible point found!"<<std::endl;
	return 0;
}

// return 1 if P is feasible and fp a point in P
// otherwise return 0 and fp has no meaning
int opt_interior(Polytope &K,
							  const int m,
							  const int n,
							  const int walk_steps,
							  const double err,
							  const double err_opt,
							  const int lw,
							  const int up,
							  const int L,
							  RNGType &rng,
							  boost::variate_generator< RNGType, boost::normal_distribution<> >
							  &get_snd_rand,
							  boost::random::uniform_real_distribution<> urdist,
							  boost::random::uniform_real_distribution<> urdist1,
							  Vector &z,
							  Vector &w){
	
	//first compute a feasible point in K
	Point fp;
  if (feasibility(K,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp)==0){
	  std::cout<<"The input polytope is not feasible!"<<std::endl;
	  return 1;
	}
	z = fp - CGAL::Origin();
	// Initialization
	Hyperplane H(fp,w);
	K.push_back(H);
	
	// create a vector V with the random points
	std::vector<Point> V;
	//initialize V !!!! This is a case without theoretical guarantees 
	for(int i=0; i<m; ++i){
		Point newv = CGAL::Origin() + z;
		hit_and_run(newv,K,n,err,rng,urdist,urdist1);
	  V.push_back(newv);
	}
	//
	int step=0;
	double len;	
	do{
		
		// compute m random points in K stored in V 
	  multipoint_random_walk(K,V,m,n,walk_steps,err,rng,get_snd_rand,urdist,urdist1);
			
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
		len = std::abs(w*newz - w*z);
		std::cout<<"step "<<step<<": "<<"z="<<z<<" "
		         "newz="<<newz<<" "
		         <<"w*z="<<w*z<<" "
		         <<"w*newz="<<w*newz<<" "
		         <<"len="<<len<<std::endl;
		
		//Update z
		z=newz;
		
		//update P with the hyperplane passing through z	
		Hyperplane H(CGAL::Origin()+z,w);
		K.pop_back();
		K.push_back(H);
		
		//check for the rest rand points which fall in new P
		std::vector<Point> newV;
		for(;vit!=V.end();++vit){
			if(Sep_Oracle(K,*vit).get_is_in())
				newV.push_back(*vit);
		}
		V=newV;
		++step;
		
	  std::cout<<"Cutting hyperplane direction="
			         <<H.orthogonal_direction()<<std::endl;
		std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
		
		if(V.empty()){
			std::cout<<"No random points left. Current OPT ="<<z<<std::endl;
			return 0;
		}
		
	}while(len>err_opt);
	
	std::cout<<"OPT = "<<z<<std::endl;
	return 0;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

int main(const int argc, const char** argv)
{ 
  double tstartall, tstopall, tstartall2, tstopall2;

  /* CONSTANTS */
  //dimension
  const size_t n=2; 
  //number of random points
  //const int m = 2*n*std::pow(std::log(n),2);
  const int m = 20*n;
  //number of walk steps
  //const int walk_steps=m*std::pow(n,3)/100;
  const int walk_steps=500; 
  //error in hit-and-run bisection of P 
  const double err=0.000001; 
  const double err_opt=0.01; 
  //bounds for the cube	
  const int lw=0, up=10000, R=up-lw;
  
  std::cout<<"m="<<m<<"\n"<<"walk_steps="<<walk_steps<<std::endl;
 
		
  /* RANDOM NUMBERS */
  // the random engine with time as a seed
  RNGType rng((double)time(NULL));
  // standard normal distribution with mean of 0 and standard deviation of 1 
	boost::normal_distribution<> rdist(0,1); 
	boost::variate_generator< RNGType, boost::normal_distribution<> > 
											get_snd_rand(rng, rdist); 
  // uniform distribution 
  boost::random::uniform_real_distribution<>(urdist); 
  boost::random::uniform_real_distribution<> urdist1(-1,1); 
  
  
  //this will store the generated random points
  std::vector<Point> V;	
	
	//this is the input polytope
	Polytope K=cube(n,lw,10);
	
	V_polytope P1, P2, Msum;
	P1.push_back(Point(-1,1));  
  P1.push_back(Point(2,1));  
  P1.push_back(Point(-1,-2));  
  
  P2.push_back(Point(1,1));  
  P2.push_back(Point(1,-1));  
  P2.push_back(Point(-1,-1));  
  P2.push_back(Point(-1,1));  
  std::cout<<!P1.empty()<<std::endl;
  std::cout<<!P2.empty()<<std::endl;
  /*
  //Transform P1, P2 to contain the origin in their interior
	Vector P1sum(n, CGAL::NULL_VECTOR);
	for(V_polytope::iterator pit=P1.begin(); pit!=P1.end(); ++pit)
		P1sum += (*pit)-CGAL::Origin();
	P1sum = P1sum/int(P1.size());
	for(V_polytope::iterator pit=P1.begin(); pit!=P1.end(); ++pit)
		*pit = CGAL::Origin() + ((*pit - CGAL::Origin()) - P1sum);
	//
	Vector P2sum(n, CGAL::NULL_VECTOR);
	for(V_polytope::iterator pit=P2.begin(); pit!=P2.end(); ++pit)
		P2sum += (*pit)-CGAL::Origin();
	P2sum = P2sum/int(P2.size());
	for(V_polytope::iterator pit=P2.begin(); pit!=P2.end(); ++pit)
		*pit = CGAL::Origin() + ((*pit - CGAL::Origin()) - P2sum);
	
  // compute mink sum using a naive algorithm
  Minkowski_sum_naive(P1,P2,Msum);
  for(V_polytope::iterator pit=Msum.begin(); pit!=Msum.end(); ++pit)
		std::cout<<*pit<<std::endl;
	*/
	
	
	//Perform optimization in the dual --> separation oracle for Minksum
	
	//Build the separation in dual 
	//query point q
	//std::cout<<"--------"<<P1sum<<" "<< P2sum<<std::endl;
	Point q(1.0/2.0,-1.0/3.0);
	//q -= (P1sum + P2sum);
	
	MinkSumPolytope P(P1,P2);
	std::cout<<Sep_Oracle(P,q).get_is_in()<<std::endl;	
	
	
	exit(1);
	
  /* OPTIMIZATION */
  //given a direction w compute a vertex v of K that maximize w*v 
  std::vector<NT> ww(n,1);
  Vector w(n,ww.begin(),ww.end());
  std::cout<<"w=";
	round_print(w);
	round_print(w/w.squared_length());
	std::cout<<w/w.squared_length()<<std::endl;
	//normalize w
	w=w/w.squared_length();
	std::cout<<"00:"<<w*Vector(0,0)<<std::endl;
	std::cout<<"10:"<<w*Vector(1,0)<<std::endl;
	std::cout<<"01:"<<w*Vector(0,1)<<std::endl;
	std::cout<<"11:"<<w*Vector(1,1)<<std::endl;
	
	
	
  const int L=30;
  
  // Interior point algorithm for optimization
  //Vector z;
	//opt_interior(K,m,n,walk_steps,err,err_opt,lw,up,L,rng,get_snd_rand,urdist,urdist1,z,w);
  
	
	
	/* Optimization with bisection
	 * 
	 */ 

  //first compute a feasible point in K (if K is non empty) 
  
  Point fp;
  optimization(K,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp,w);
  std::cout<<"OPT="<<fp<<std::endl;
  /*
  if (feasibility(K,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp)==0){
	  std::cout<<"The input polytope is not feasible!"<<std::endl;
	  return 1;
	}
	
	//then compute a point outside K along the line (fp,w)
  Point pout=fp+100*w;
  Point pin=fp;
  
  
  std::cout<<"Start point: ";
  round_print(pout);
  Vector aug(w);
  while(Sep_Oracle(K,pout).get_is_in() == true){
    aug*=2;
    pout+=aug;
    std::cout<<"Next point: ";
    round_print(pout);
  }
  
  //find a hyperplane that is not feasible
  bool feasible=true;
  do{
	  Hyperplane H(pout,w);
	  std::cout<<std::endl<<"CHECKING FEASIBILITY IN :"<<pout<<std::endl;
		K.push_back(H);
		if(feasibility(K,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp) == 1){
			aug*=2;
      pout+=aug;
      std::cout<<"Outside point but feasible hyperplane: ";
    }
		else
			feasible=false;
    K.pop_back();
  }while(feasible);
  std::cout<<"NON feasible hyperplane found. pout= ";
  round_print(pout);
  
  
  //binary search for optimization
  double len;
  Point pmid;
  do{
		pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
		Hyperplane H(pmid,w);
		K.push_back(H);
		std::cout<<"pmid,pin,pout,w"<<std::endl;
		round_print(pmid);
		round_print(pin);round_print(pout);round_print(w);
		
		if(feasibility(K,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp) == 1)
			pin=pmid;
		else
			pout=pmid;
		K.pop_back();
		len=std::abs((pin-CGAL::Origin())*w - (pout-CGAL::Origin())*w);
		std::cout<<"len="<<len<<std::endl;
		std::cout<<"fp=";round_print(fp);
	}while(len > err_opt);
	std::cout<<"fp=";
	round_print(fp);
	std::cout<<"w=";
	round_print(w);
	*/
	
  return 0;
}
