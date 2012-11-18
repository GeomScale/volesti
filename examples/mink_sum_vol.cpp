// Copyright 2012-2013 National and Kapodistrian University of Athens, Greece.
//
// This file is part of RandGeom.
//
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

#include <geom_rand_walks.h>

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

int main(const int argc, const char** argv)
{ 
	// VARS
	int n;
	double wl_c, e;
	//parse command line input
	if(argc==4){
		std::cout << argv[3]<<std::endl;
	  //dimension
		n = atoi(argv[1]);
		//constants
		e = atof(argv[2]);
		wl_c = atof(argv[3]);
	}else{
		std::cout<<"Wrong num of args"<<std::endl;
		exit(1);
	}
	
  double tstart, tstop;

  /* CONSTANTS */
  //int n=4; 
  //number of random points
  //const int m = 2*n*std::pow(std::log(n),2);
  const int m = 20*n;
  //number of walk steps
  //const int walk_steps=m*std::pow(n,3)/100;
  const int walk_steps=20*n; 
  //
  const int L=30;
  //error in hit-and-run bisection of P 
  const double err=0.000001; 
  const double err_opt=0.01; 
  //bounds for the cube	
  const int lw=0, up=10000, R=up-lw;
  
  //std::cout<<"m="<<m<<"\n"<<"walk_steps="<<walk_steps<<std::endl;
 
		
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
	


  /*
  Hyperplane H(n,fp.cartesian_begin(),fp.cartesian_end(),1);
  
  std::cout<<"which side(P)="<<H.has_on_positive_side(CGAL::Origin()+w)<<std::endl;
  std::cout<<"which side(N)="<<H.has_on_negative_side(CGAL::Origin()+w)<<std::endl;
  
  Point q(1.0/2.0,-1.0/3.0);
  std::cout<<"Test"<<std::endl;
	std::cout<<"is in:"<<Sep_Oracle(Msum,q).get_is_in()<<std::endl;	
	std::cout<<"H sep:"<<Sep_Oracle(Msum,q).get_H_sep()<<std::endl;	
  Hyperplane Htest = Sep_Oracle(Msum,q).get_H_sep();
  std::cout<<"H dim:"<<Htest.dimension()<<std::endl;	
  for(Hyperplane::Coefficient_const_iterator Hit=Htest.coefficients_begin();
								Hit!=Htest.coefficients_end(); ++Hit)
		std::cout<<*Hit<<" ";
	std::cout<<std::endl;
	*/
	
	
  
  
  // VOLUME of MINKOWSKI SUM 
  
  V_polytope P1, P2;
	P1.push_back(Point(-1,1));  
  P1.push_back(Point(2,1));  
  P1.push_back(Point(-1,-2));  
  
  P2.push_back(Point(1,1));  
  P2.push_back(Point(1,-1));  
  P2.push_back(Point(-1,-1));  
  P2.push_back(Point(-1,1));  
  std::cout<<!P1.empty()<<std::endl;
  std::cout<<!P2.empty()<<std::endl;
  
  MinkSumPolytope Msum(P1,P2);
  
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
	Vector q2(2,-2);
	//q -= (P1sum + P2sum);
	
	
	Point fp;
  optimization(Msum,m,n,walk_steps,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp,q2);
  std::cout<<"OPT="<<fp<<std::endl;
  
  Hyperplane H(n,fp.cartesian_begin(),fp.cartesian_end(),1);
  
  //std::cout<<"prod="<<(fp-CGAL::Origin())*q2<<std::endl;
  if(double((fp-CGAL::Origin())*q2) <= double(1.0))
    std::cout<<"result=IN"<<std::endl;
  else
    std::cout<<"result=OUT"<<std::endl;
  
  std::cout<<"which side(P)="<<H.has_on_positive_side(CGAL::Origin()+q2)<<std::endl;
  std::cout<<"which side(N)="<<H.has_on_negative_side(CGAL::Origin()+q2)<<std::endl;
  
  Point q(1.0/2.0,-1.0/3.0);
  std::cout<<"Test"<<std::endl;
	std::cout<<"is in:"<<Sep_Oracle(Msum,q).get_is_in()<<std::endl;	
	std::cout<<"H sep:"<<Sep_Oracle(Msum,q).get_H_sep()<<std::endl;	
  Hyperplane Htest = Sep_Oracle(Msum,q).get_H_sep();
  std::cout<<"H dim:"<<Htest.dimension()<<std::endl;	
  for(Hyperplane::Coefficient_const_iterator Hit=Htest.coefficients_begin();
								Hit!=Htest.coefficients_end(); ++Hit)
		std::cout<<*Hit<<" ";
	std::cout<<std::endl;

  

	
  return 0;
}


