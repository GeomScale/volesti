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

// ball type
struct Ball{
	public:
	  Ball(Point c, NT R) : _c(c),	 _R(R) {}
	   
	  Point center(){
			return _c;
		}
		NT squared_radius(){
			return _R;
		}
		NT radius(){
			return std::sqrt(_R);
		}
		bool is_in(Point p){
			return ((p-CGAL::Origin()) - (_c-CGAL::Origin())).squared_length() <= _R;
		}
	
	private:
	  Point  _c; //center
	  NT     _R; //squared radius
};


// define different kind of polytopes
typedef std::vector<Hyperplane>           H_polytope;
typedef H_polytope 								        Polytope;
typedef std::vector<Point>						    V_polytope;

typedef std::pair<V_polytope,V_polytope> 	    MinkSumPolytope;
typedef std::pair<MinkSumPolytope,bool> 	    MinkSumPolytopeDual;

template <class T>
class BallIntersectPolytope {
  public:
    BallIntersectPolytope(T &P, Ball &B) : _P(P), _B(B) {};
    
    T first() { return _P; }
    Ball second() { return _B; }
    
  private:  
    T    _P;
    Ball _B;
};


// define random generators
typedef boost::mt19937 RNGType; ///< mersenne twister generator
//typedef boost::lagged_fibonacci607 RNGType;

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
				  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
	  boost::variate_generator< RNGType, boost::normal_distribution<> >
	  &get_snd_rand;
	  boost::random::uniform_real_distribution<> urdist;
	  boost::random::uniform_real_distribution<> urdist1;
};

// define extreme points
typedef CGAL::Extreme_points_traits_d<Point>   EP_Traits_d;


//function to print rounding to double coordinates 
template <class T>
void round_print(T p) { 
  for(typename T::Cartesian_const_iterator cit=p.cartesian_begin(); 
      cit!=p.cartesian_end(); ++cit)
	  std::cout<<CGAL::to_double(*cit)<<" "; 
  std::cout<<std::endl;
}

// Naive algorithm for Mink sum
typedef std::vector<V_polytope>              Vpolys;

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
Polytope cube(int n, NT lw, NT up){	
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
		//std::cout<<apex[i]<<" ";
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

// contruct a n-ball of radius r centered in the origin  
/*
Ball ball(int n, const NT r){
  
  std::vector<Point> P_ball;
  for(int i=0; i<n; ++i){
		std::vector<NT> coords;
		for(int j=0; j<n; ++j){
			if(i==j) 
				coords.push_back(r);
			else coords.push_back(NT(0));
		}
		P_ball.push_back(Point(n,coords.begin(),coords.end()));
	}
	std::vector<NT> extra_coords(n,NT(0));
	extra_coords[0]=NT(-1*r);
	P_ball.push_back(Point(n,extra_coords.begin(),extra_coords.end()));
  Ball B(n,P_ball.begin(),P_ball.end());
	return B;
}
*/

//template <typename T> struct Oracle{
//  sep Sep_Oracle(T &P, Point v);
//};

template <class T>
int optimization(T &KK,vars var,Point &fp,Vector &w);
template <class T>
int opt_interior(T &K,vars &var,Point &opt,Vector &w);

// GENERIC ORACLE DESCRIPTION
// template parameter specialization
template <typename T> 
sep Sep_Oracle(T&, Point, vars&);


//function that implements the separation oracle 
template<> sep Sep_Oracle<Polytope>(Polytope &P, 
                                    Point v,
                                    vars &var)
{
	typename Polytope::iterator Hit=P.begin(); 
	while (Hit!=P.end()){
		if (Hit->has_on_negative_side(v))
			return sep(false,*Hit);
		++Hit;
	}
	return sep(true);	
}

// BallIntersectPolytope separation oracle
template <typename T> sep Sep_Oracle(BallIntersectPolytope<T> &P, 
                                     Point v,
                                     vars &var)
{
	T P1 = P.first(); 
	Ball B = P.second();
	
	if (B.is_in(v)){
		return Sep_Oracle(P1,v,var);
	}
	// the problem here is that it is out but without separating hyperplane
	// so this is a membership oracle! not separation
	// TODO: fix it!
	return sep(false);
}

// Minkowski sum Separation 
// this is Optimization in the dual
template<> sep Sep_Oracle<MinkSumPolytope>(MinkSumPolytope &P, 
																				   Point v,
																				   vars &var)
{	
	MinkSumPolytopeDual Pdual(P,true);
	Point fp;
	Vector w=v-CGAL::Origin();
  double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
  //std::cout<<"SEP_ORACLE #randpoints="<<var.m<<std::endl;
  optimization(Pdual,var,fp,w);
  //std::cout<<"O"<<(double)clock()/(double)CLOCKS_PER_SEC - tstart<<" "<<std::flush;
  //opt_interior(Pdual,var,fp,w);
  //std::cout<<"OPT="<<fp<<std::endl;
  if ((fp-CGAL::Origin()) * (v-CGAL::Origin()) <= NT(1.0))
    return sep(true); // is in
  else { // is out
    Hyperplane H_sep(fp,(v-CGAL::Origin()).direction());
    return sep(false,H_sep);
  }
}	


// DUAL Minkowski sum Separation 
// this is Optimization in the primal
template<> sep Sep_Oracle<MinkSumPolytopeDual>(MinkSumPolytopeDual &Pdual, 
																				       Point v,
																				       vars &var)
{	
	MinkSumPolytope P = Pdual.first;
	//tranform query point v to a vector q
	Vector q=v-CGAL::Origin();
	V_polytope P1=P.first;
	V_polytope P2=P.second;
	NT max = q * (*(P1.begin())-CGAL::Origin());
	Point max_p1 = *(P1.begin());
	for(V_polytope::iterator pit=P1.begin(); pit!=P1.end(); ++pit){
		double innerp = q * (*pit-CGAL::Origin());
		//std::cout<<*pit<<" "<<q<<" "<<innerp<<std::endl;
		if(innerp > max){
			max = innerp;
			max_p1 = *pit;
		}
	}
	//std::cout<<max_p1<<std::endl;
	max = q * (*(P2.begin())-CGAL::Origin());
	Point max_p2 = *(P2.begin());
	for(V_polytope::iterator pit=P2.begin(); pit!=P2.end(); ++pit){
		double innerp = q * (*pit-CGAL::Origin());
		if(innerp > max){
			max = innerp;
			max_p2 = *pit;
		}
	}
	//std::cout<<max_p2<<std::endl;
	Vector max_psum=(max_p1-CGAL::Origin())+(max_p2-CGAL::Origin());
	//+ q*(-P1sum -P2sum)
	if(max_psum*q <= 1)
		return sep(true);	
	else{
		//construct the dual hyperplane of max_psum
	  //std::vector<NT> vcoord;
	  //for(Vector::Cartesian_const_iterator cit=max_psum.cartesian_begin(); 
    //    cit!=max_psum.cartesian_end(); ++cit){
		//  vcoord.push_back(*cit);
		//}
	  Hyperplane H_sep(v,-max_psum.direction());
	  return sep(false,H_sep);
	}
	//if(max_psum*q == 1)
	//	std::cout<<"sharp!"<<std::endl;
}	
 
 
// function to find intersection of a line and a polytope 
template <class T>
Vector line_intersect(Point pin, 
                      Vector l, 
                      T &P, 
											vars &var){
	
	double err = var.err;											
  Vector vin=pin-CGAL::Origin();
  //first compute a point outside P along the line
  Point pout=pin;
  //std::cout<<"Starting inside point: "<<vin<<std::endl;
	
	Vector aug(l);
  while(Sep_Oracle(P,pout,var).get_is_in() == true){
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
		if(Sep_Oracle(P,CGAL::Origin()+vmid,var).get_is_in() == false)
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

//version 2
template <class T>
Vector line_intersect(Point pin, 
                      Vector l, 
                      T &P, 
											vars &var,
											vars &var2){
	
	double err = var.err;											
  Vector vin=pin-CGAL::Origin();
  //first compute a point outside P along the line
  Point pout=pin;
  //std::cout<<"Starting inside point: "<<vin<<std::endl;
	
	Vector aug(l);
  while(Sep_Oracle(P,pout,var2).get_is_in() == true){
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
		if(Sep_Oracle(P,CGAL::Origin()+vmid,var2).get_is_in() == false)
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

/* Hit and run Random Walk */
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
	
	std::vector<NT> v;
	for(int i=0; i<n; ++i)
		v.push_back(urdist1(rng));
	Vector l(n,v.begin(),v.end());
	Vector b1 = line_intersect(p,l,P,var);
	Vector b2 = line_intersect(p,-l,P,var);
	//std::cout<<"b1="<<b1<<"b2="<<b2<<std::endl;
	double lambda = urdist(rng);
	p = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
	return 1;
}

//version 2
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
	
	std::vector<NT> v;
	for(int i=0; i<n; ++i)
		v.push_back(urdist1(rng));
	Vector l(n,v.begin(),v.end());
	Vector b1 = line_intersect(p,l,P,var,var2);
	Vector b2 = line_intersect(p,-l,P,var,var2);
	//std::cout<<"b1="<<b1<<"b2="<<b2<<std::endl;
	double lambda = urdist(rng);
	p = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
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
	 boost::variate_generator< RNGType, boost::normal_distribution<> >
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
			Vector b1=line_intersect(v,l,P,var);	
			Vector b2=line_intersect(v,-l,P,var);
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
	const int m = var.m;
  int n = var.n;
  const int walk_steps = var.walk_steps;
  const double err = var.err;
  const int lw = var.lw;
  double up = var.up;
  const int L = var.L;
  RNGType &rng = var.rng;
  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
			if (print) std::cout<<"Feasible point found! "<<z<<std::endl;
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
	if (print) std::cout<<"No feasible point found!"<<std::endl;
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
  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
    if (print2) std::cout<<"Step:"<<step<<" Current fp= "<<fp
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
  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
		if (print2) std::cout<<"Step:"<<step<<" Current centroid= "<<z
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
	bool print2 = true;
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
  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
	int m = var.m;
	int n = var.n;
	int walk_steps = var.walk_steps;
	const double err = var.err;
	const double err_opt = var.err_opt;
	const int lw = var.lw;
	double up = var.up;
	const int L = var.L;
  RNGType &rng = var.rng;
  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
		
		if (print) 
		  std::cout<<"Step:"<<step<<" Current centroid= "<<z
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
	
	bool print = false;
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
  for(int i=0; i<rnum; ++i){ //generate rnum rand points 
		//start with a u.d.r point in the smallest ball 
		//radius=1, center=Origin()
		std::vector<NT> coords(n,0);
		Point p(n,coords.begin(),coords.end());
		BallPoly PBold(P,balls[0]);
		hit_and_run(p,PBold,var,var2);
		//std::cout<<p<<std::endl;
		//std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
		//std::cout<<balls[0].is_in(p)<<std::endl;
		
		//exit(0);
		std::vector<Ball>::iterator bit=balls.begin();
		std::vector<int>::iterator prod_it=telescopic_prod.begin();
		++bit; 
		if (print) std::cout<<"\n\ngenerate random point..."<<i<<"/"<<rnum<<" ";
		const NT pi = boost::math::constants::pi<NT>();
		NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0); 
		for(std::vector<int>::iterator prod_it=telescopic_prod.begin(); 
		    prod_it!=telescopic_prod.end(); ++prod_it){
			vol *= NT(i+1)/NT(*prod_it);
		}
		if (print) std::cout<<"current vol estimation= "<<vol<<std::endl;
	  if (print) std::cout<<"walklen="<<walk_len<<std::endl;
	  if (print) std::cout<<"rnum="<<rnum<<std::endl;
		for(; bit!=balls.end(); ++bit, ++prod_it){
			// generate a random point in bit intersection with P 
			BallPoly PB(P,*bit);
			
			for(int j=0; j<walk_len; ++j){
			  hit_and_run(p,PB,var,var2);
				//std::cout<<"h-n-r:"<<p<<std::endl;
			}
			if (Sep_Oracle(PBold,p,var2).get_is_in()){
				//std::cout<<p<<" IN ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
			  ++(*prod_it);
			}else{
			  ;
			  //std::cout<<p<<":"<<(p-CGAL::Origin()).squared_length()
			  //<<" OUT ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
			}
			PBold=PB;
		}
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
  boost::variate_generator< RNGType, boost::normal_distribution<> >
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
	os << "print double($f);\n";
	os << "print ' ';\n";
	os << "print tv_interval($t0,[gettimeofday]);\n";
	os << "print ' ';\n";
	os << std::endl;
}

