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


typedef CGAL::Cartesian_d<NT> 	      Kernel; 
//typedef CGAL::Triangulation<Kernel> T;
typedef Kernel::Point_d								Point;
typedef Kernel::Vector_d							Vector;
typedef Kernel::Line_d								Line;
typedef Kernel::Hyperplane_d					Hyperplane;
typedef Kernel::Direction_d						Direction;
typedef std::vector<Hyperplane>       H_polytope;
typedef H_polytope 								    Polytope;

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
													 boost::random::uniform_real_distribution<> &urdist){
													
	//generate more points (using points in V) in order to have m in total
	std::vector<Point> U;
	for(int mk=0; mk<m-V.size(); ++mk){
		// Compute a point as a random uniform convex combination of V 
		std::vector<double> a;
		double suma=0;
		for(int ai=0; ai<V.size(); ++ai){
			a.push_back(urdist(rng));
			suma+=a[a.size()-1];
		}		
		std::vector<Point>::iterator Vit=V.begin();
		Vector p(n,CGAL::NULL_VECTOR);
		for(std::vector<double>::iterator ait=a.begin(); ait!=a.end(); ++ait){
		  p+=NT(*ait)/NT(suma)*((*Vit)-(CGAL::Origin()));
		  ++Vit;
		}
		U.push_back(CGAL::Origin()+p);
	}
	//append U to V
	V.insert(V.end(),U.begin(),U.end());											 
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
			//round_print(v);
			*vit=v;
	  }
	}
}

// return 1 if P is feasible and fp a point in P
// otherwise return 0 and fp has no meaning
int feasibility(Polytope &KK,
                std::vector<Point> &V,
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
							  boost::random::uniform_real_distribution<> &urdist,
							  Point &fp){
	
	//this is the large cube contains the polytope
	Polytope P=cube(n,lw,up);								
  int step=0;
  while(step < 2*n*L){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,m,n,walk_steps,err,rng,get_snd_rand,urdist);
		
	  //compute the average using the half of the random points
		Vector z(n,CGAL::NULL_VECTOR);
		int i=0;
		std::vector<Point>::iterator vit=V.begin();
		for(; i<m/2; ++i,++vit){
			CGAL:assert(vit!=V.end());
			z = z + (*vit - CGAL::Origin());
		}
		z=z/m;	
		std::cout<<"step "<<step<<": "<<"z=";
		round_print(z);
		
		sep sep_result = Sep_Oracle(KK,CGAL::Origin()+z);
		if(sep_result.get_is_in()){
			std::cout<<"Feasible point found!"<<std::endl;
			fp = CGAL::Origin() + z;
			return 1;
		}
		else {
			//update P with the hyperplane passing through z
			Hyperplane H(CGAL::Origin()+z,sep_result.get_H_sep().orthogonal_direction());
			P.push_back(H);
			//check for the rest rand points which fall in new P
			std::vector<Point> newV;
			for(;vit!=V.end();++vit){
				if(Sep_Oracle(P,*vit).get_is_in())
					newV.push_back(*vit);
			}
			V=newV;
			++step;
		}
	}
	std::cout<<"No feasible point found!"<<std::endl;
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
  const size_t n=3; 
  //number of random points
  const int m=2*30;
  //number of walk steps
  const int walk_steps=100;
  //error in hit-and-run bisection of P
  const double err=0.000001;
  const double err_opt=0.000001;  
  //bounds for the cube
  const int lw=0, up=10000, R=up-lw;
  
  /* INITIALIZE POINTS IN CUBE*/ 
  CGAL::Random CGALrng;
  std::vector<Point> V;
  for(size_t i=0; i<m; ++i){
		std::vector<NT> t;
		for(size_t j=0; j<n; ++j)
			t.push_back(NT(CGALrng.get_int(lw,up)));
		Point v(n,t.begin(),t.end());
		V.push_back(v);
		//std::cout<<v<<std::endl;
	}
		
	
	//this is the polytope
	Polytope K=cube(n,lw,10);
	
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
  
  /* OPTIMIZATION */
  //given a direction w compute a vertex v of K that maximize w*v 
  Vector w(*(V.begin())-CGAL::Origin());
  std::cout<<"w=";
	round_print(w);
	
  const int L=10;
  //first compute a feasible point in K (if K is non empty) 
  Point fp;
  feasibility(K,V,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,fp);	  
	
	//then compute a point outside K along the line (fp,w)
  Point pout=fp;
  Point pin=fp;
  
  Vector aug(w);
  while(Sep_Oracle(K,pout).get_is_in() == true){
    aug*=2;
    pout+=aug;
    //std::cout<<"Outside point: ";
    //round_print(pout);
  }
  
  //binary search for optimization
  double len;
  Point pmid;
  do{
		pmid=CGAL::Origin()+(((pin-CGAL::Origin())+(pout-CGAL::Origin()))/2);
		Hyperplane H(pmid,w);
		K.push_back(H);
		round_print(pmid);
		
		if(feasibility(K,V,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,fp) == 1)
			pin=pmid;
		else
			pout=pmid;
		K.pop_back();
		len=CGAL::to_double(((pin-CGAL::Origin())+(pout-CGAL::Origin())).squared_length());
		std::cout<<"len="<<len<<std::endl;
	}while(len > err_opt);
	std::cout<<"fp=";
	round_print(fp);
	std::cout<<"w=";
	round_print(w);
	
	
  return 0;
}
