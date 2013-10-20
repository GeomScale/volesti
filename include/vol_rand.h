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
#include <ballintersectpolytope.h>
//#include <opt_rand.h>
//#include <oracles.h>
#include <random_samplers.h>
#include <misc.h>

/////////////////////////////////////////////////////////
// VOLUME
// randomized approximate volume computation 
/*************************************************
/* VOLUME with random DIRECTIONS hit and run     */
// We assume that the polytope P is properly sandwitched
// The sandwitching:
// r is the radius of the smallest ball
// d is the radius of the largest
template <class T>
NT volume0(T &P,
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


/*************************************************
/* VOLUME with random COORDINATES hit and run    */
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
  boost::random::uniform_int_distribution<> uidist(0,n-1);
  //boost::random::uniform_real_distribution<> urdist1 = var.urdist1;
	
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
		CGAL::Random_points_in_ball_d<Point> gen (n, NT(1.0));
		p = *gen;
		//std::cout<<p<<std::endl;
		//std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
		//std::cout<<balls[0].is_in(p)<<std::endl;
		
		std::vector<Ball>::iterator bit=balls.begin();
		std::vector<int>::iterator prod_it=telescopic_prod.begin();
		++bit; 
		for(; bit!=balls.end(); ++bit, ++prod_it){
			// generate a random point in bit intersection with P 
			BallPoly PB(P,*bit);
			
			std::vector<NT> lamdas(P.size(),NT(0));
			int rand_coord = uidist(rng);
			double kapa = urdist(rng);
			Point p_prev=p;
  		hit_and_run_coord_update(p,p_prev,PB,rand_coord,rand_coord,kapa,lamdas,var,var2,true);
							
			for(int j=0; j<walk_len; ++j){
			  int rand_coord_prev = rand_coord;
			  rand_coord = uidist(rng);
			  kapa = urdist(rng);
				hit_and_run_coord_update(p,p_prev,PB,rand_coord,rand_coord_prev,kapa,lamdas,var,var2,false);
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
	  //if (print) std::cout<<"rnum="<<rnum<<std::endl;
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

