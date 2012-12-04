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

// Approximating the volume of the Minkowski sum of convex 
// polytopes or bodies. 
// The user should provide the appropriate membership 
// oracles for the summands.

int main(const int argc, const char** argv)
{ 
	// VARS
	int n;
	double wl_c, e;
	//parse command line input
	if(argc==4){
		//std::cout << argv[3]<<std::endl;
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
  int m = 20*n;
  //number of walk steps
  //const int walk_steps=m*std::pow(n,3)/100;
  const int walk_steps=20*n; 
  //
  const int L=30;
  //error in hit-and-run bisection of P 
  const double err=0.001; 
  const double err_opt=0.0001; 
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
  
  
  // VOLUME of MINKOWSKI SUM 
 

  // for experiments
  //for (n=2; n<8; ++n){
  
  // Input polytopes
  //
  V_polytope P1, P2;
	//std::cout<<"cube"<<std::endl;
	P1=Vcube(n,-1,1);
	//std::cout<<"cross"<<std::endl;
	P2=Vcross(n,-1,1);

	/* 2d example
	P1.push_back(Point(-1,1));  
  P1.push_back(Point(2,1));  
  P1.push_back(Point(-1,-2));  
  
  P2.push_back(Point(1,1));  
  P2.push_back(Point(1,-1));  
  P2.push_back(Point(-1,-1));  
  P2.push_back(Point(-1,1));  
  std::cout<<!P1.empty()<<std::endl;
  std::cout<<!P2.empty()<<std::endl;
  */

  MinkSumPolytope Msum(P1,P2);
  
  /* AFFINE TRANFORMATION OF THE POLYTOPES TO CONTAIN ORIGIN 
   
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
	
	//volume vars
  // (#rand_points, n, #walk_steps, ...)
  int rnum = std::pow(e,-2) * 400 * n * std::log(n);
  int walk_len =  wl_c * std::pow(n,4);
  //rnum=e;
  //walk_len=wl_c;

  vars var1(rnum,n,walk_len,err,err_opt,lw,1,L,rng,get_snd_rand,urdist,urdist1);
  
  //opt vars
  // (#rand_points, n, #walk_steps, ...)
  int rnum_opt =  5 * n * std::pow(std::log(n),2);
  int walk_len_opt =  0.001 * std::pow(n,4);

  vars var2(rnum_opt,n,walk_len_opt,err,err_opt,lw,1,L,rng,get_snd_rand,urdist,urdist1);
  
  //RUN THE EXPERIMENTS
  int num_of_exp=10;
  for(int i=0; i<num_of_exp; ++i){
	  tstart = (double)clock()/(double)CLOCKS_PER_SEC;
	  double v1 = volume1(Msum,var1,var2,1,std::sqrt(5.0));
	  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
	  //double v2 = volume2(P,n,rnum,walk_len,err,rng,get_snd_rand,urdist,urdist1);
	  /*
	  double exactvol = std::pow(2,n);
	  
	  std::cout<<rnum<<"\n\n\nALGORITHM 1\n-----------\nvolume = "
	           //<<(1-e)*exactvol<<" < "<<v1<<" < "<<(1+e)*exactvol
	           <<v1<<std::endl;
		//std::cout<<"exact volume = "<<exactvol<<std::endl;
		std::cout<<"# walk steps = "<<walk_len<<std::endl;
		std::cout<<"# rand points = "<<rnum<<std::endl;
		std::cout<<"time = "<<tstop-tstart<<std::endl;
		*/
		std::cout<<n<<" "
			       <<rnum<<" "
			       <<walk_len<<" "
			       <<v1<<" "
			       <<tstop-tstart<<std::endl;
	}	
  //}
	
	//EXACT ALGORITHM
	//compute the Minkowski sum and then the volume
	// compute mink sum using a naive algorithm
  /*
  V_polytope P;
  Minkowski_sum_naive(P1,P2,P);
  //for(V_polytope::iterator pit=P.begin(); pit!=P.end(); ++pit)
		//std::cout<<*pit<<std::endl;
	
	std::ofstream polymakefile;
	polymakefile.open("volume.polymake");
	print_polymake_volfile(P,polymakefile);
	system ("polymake volume.polymake");
	std::cout<<std::endl;
	*/
	
  return 0;
}


