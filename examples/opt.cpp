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
	std::cout<<"--------------------------"<<std::endl;
	std::cout<<"Random points before walk"<<std::endl;
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit){
		Point v=*vit;
		hit_and_run(v,P,n,err,rng,urdist,urdist1);
		//std::cout<<*vit<<"---->"<<v<<std::endl;
	}
	std::cout<<"WALKING......"<<std::endl;											 
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
	std::cout<<"Random points after walk"<<std::endl;
	for(std::vector<Point>::iterator vit=V.begin(); vit!=V.end(); ++vit)
		std::cout<<*vit<<std::endl;											 
	std::cout<<"--------------------------"<<std::endl;
	for(Polytope::iterator polyit=P.begin(); polyit!=P.end(); ++polyit)
		std::cout<<*polyit<<std::endl;
	
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
							  boost::random::uniform_real_distribution<> urdist,
							  boost::random::uniform_real_distribution<> urdist1,
							  Point &fp){
	
	//this is the large cube contains the polytope
	Polytope P=cube(n,lw,up);								
  int step=0;
  while(step < 2*n*L){
	  // compute m random points in P stored in V 
	  multipoint_random_walk(P,V,m,n,walk_steps,err,rng,get_snd_rand,urdist,urdist1);
		
	  //compute the average using the half of the random points
		Vector z(n,CGAL::NULL_VECTOR);
		int i=0;
		std::vector<Point>::iterator vit=V.begin();
		std::cout<<"RANDOM POINTS"<<std::endl;
		for(; i<m/2; ++i,++vit){
			CGAL:assert(vit!=V.end());
			z = z + (*vit - CGAL::Origin());
			std::cout<<*vit<<std::endl;
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
			std::cout<<"Cutting hyperplane direction="<<H.orthogonal_direction()<<std::endl;
			std::cout<<"Number of random points in new P="<<newV.size()<<"/"<<m/2<<std::endl;
			
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
  const size_t n=2; 
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
	
	Point ptest(NT(9),NT(2)); 
  Vector wtest(NT(1),NT(0)); 
  Hyperplane Htest(ptest,wtest);
  K.push_back(Htest);
  std::cout<<"hyperplane="<<Htest<<std::endl;
	
	Point ptest2(NT(9.5),NT(2)); 
  
	std::cout<<"Intersection: "<< line_intersect(ptest2,wtest,K,err)
	         <<","<<line_intersect(ptest2,-wtest,K,err)<<std::endl;
	
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
  boost::random::uniform_real_distribution<> urdist1(-1,1); 
  
  /* OPTIMIZATION */
  //given a direction w compute a vertex v of K that maximize w*v 
  Vector w(*(V.begin())-CGAL::Origin());
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
  //first compute a feasible point in K (if K is non empty) 
  Point fp;
  if (feasibility(K,V,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp)==0){
	  std::cout<<"The input polytope is not feasible!"<<std::endl;
	  return 1;
	}
	
	//then compute a point outside K along the line (fp,w)
  Point pout=fp;
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
		K.push_back(H);
		if(feasibility(K,V,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp) == 1){
			aug*=2;
      pout+=aug;
      std::cout<<"Outside point but feasible hyperplane: ";
    }
		else
			feasible=false;
    round_print(pout);
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
		std::cout<<"pmid,pin,pout"<<std::endl;
		round_print(pmid);
		round_print(pin);round_print(pout);
		
		if(feasibility(K,V,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp) == 1)
			pin=pmid;
		else
			pout=pmid;
		K.pop_back();
		len=CGAL::to_double(((pin-CGAL::Origin())+(pout-CGAL::Origin())).squared_length());
		std::cout<<"len="<<len<<std::endl;
		std::cout<<"fp=";round_print(fp);
	}while(len > err_opt);
	std::cout<<"fp=";
	round_print(fp);
	std::cout<<"w=";
	round_print(w);
	
	
  return 0;
}
