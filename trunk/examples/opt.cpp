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

//#include <gmpxx.h>
//typedef mpq_class NT;
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
//typedef CGAL::Gmpq                NT;
typedef double                NT;
//typedef CGAL::Gmpz                NT;


typedef CGAL::Cartesian_d<NT> 	  K; 
//typedef CGAL::Triangulation<K> 	T;
typedef K::Point_d								Point;
typedef K::Vector_d								Vector;
typedef K::Line_d								  Line;
typedef K::Hyperplane_d						Hyperplane;
typedef K::Direction_d						Direction;
typedef std::vector<Hyperplane>   H_polytope;
typedef H_polytope 								Polytope;

typedef boost::mt19937 RNGType; ///< mersenne twister generator

//function to print rounding to double coordinates 
template <class T>
void round_print(T p) { 
  for(typename T::Cartesian_const_iterator cit=p.cartesian_begin(); 
      cit!=p.cartesian_end(); ++cit)
	  std::cout<<CGAL::to_double(*cit)<<" "; 
  std::cout<<std::endl;
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

//function that implements the separation oracle 
template<typename Polytope> sep Sep_Oracle(Polytope P, Point v)
{
    typename Polytope::iterator Hit=P.begin(); 
    while (Hit!=P.end()){
      if (Hit->has_on_negative_side(v))
				return sep(false,*Hit);
      ++Hit;
    }
		return sep(true);	
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

/*-------- MULTIPOINT RANDOM WALK -------*/
int multipoint_random_walk(Polytope &P,
													 std::vector<Point> &V,
													 const int n,
													 const int walk_steps,
													 const double err,
													 RNGType &rng,
													 boost::variate_generator< RNGType, boost::normal_distribution<> >
													 &get_snd_rand,
													 boost::random::uniform_real_distribution<> &urdist){
													
													 
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
			
			// Move the point to a random (uniform) point in P along the constructed line 
			double lambda = urdist(rng);		
			v = CGAL::Origin() + (NT(lambda)*b1 + (NT(1-lambda)*b2));
			//std::cout<<v<<std::endl;
			round_print(v);
			*vit=v;
	  }
	}
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
  const size_t n=3; 
  //number of random points
  const int m=100;
  //number of walk steps
  const int walk_steps=1000;
  //error in hit-and-run bisection of P
  const double err=0.000001;
  //bounds for the cube
  const int lw=0, up=1000;
  
  /* INITIALIZE POINTS */ 
  CGAL::Random CGALrng;
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(CGALrng.get_int(lw,up)));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		std::cout<<v<<std::endl;
	}
		
	Polytope P=cube(n,lw,up);
	Polytope KP=cube(n,lw,10);
	
	//compute the average
	Vector z0(n,CGAL::NULL_VECTOR);
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit){
		z0 = z0 + (*vit - CGAL::Origin());
	}	
	z0=z0/m;
	std::cout<<"z=";
	round_print(z0);

  // RANDOM NUMBERS
  // the random engine with time as a seed
  RNGType rng((double)time(NULL));
  // standard normal distribution with mean of 0 and standard deviation of 1 
	boost::normal_distribution<> rdist(0,1); 
	boost::variate_generator< RNGType, boost::normal_distribution<> >
											get_snd_rand(rng, rdist); 
  // uniform distribution
  boost::random::uniform_real_distribution<>(urdist); 
  
  // compute m random points in P stored in V 
  multipoint_random_walk(P,V,n,walk_steps,err,rng,get_snd_rand,urdist);
	
  //compute the average
	Vector z(n,CGAL::NULL_VECTOR);
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit){
		z = z + (*vit - CGAL::Origin());
	}
	z=z/m;	
	std::cout<<"z=";
	round_print(z);
	  
	std::cout<<"z0=";
	round_print(z0);
	
  return 0;
}
