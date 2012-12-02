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

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// Approximating the volume of a convex polytope or body 
// can also be used for integration of concave functions.
// The user should provide the appropriate membership 
// oracles.

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
	
	//timings
  double tstart, tstop;


  /* CONSTANTS */
  //error in hit-and-run bisection of P 
  const double err=0.000001; 
  const double err_opt=0.01; 
  //bounds for the cube	
  const int lw=0, up=10000, R=up-lw;
  
 
		
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


	//for(n=2; n<13; n+=2){
	//std::cout<<"n="<<n<<"\n"<<std::endl;
	/* VOLUME */
  
  /* CUBE 
	Polytope P = cube(n,-1,1);
  //sandwitch
  std::vector<NT> coords_apex(n,1);
	Vector p_apex(n,coords_apex.begin(),coords_apex.end());
  NT r=1, d=std::sqrt(p_apex.squared_length());
  /**/
  
  /*  Cross-polytope */
  Polytope P=cross(n,-1,n);
  std::vector<NT> apex;
	for(int j=0; j<n; ++j){
		if(j==0) apex.push_back(NT(n));
		else apex.push_back(NT(0));
  }
  Vector p_apex(n,apex.begin(),apex.end());
  NT r=1, d=std::sqrt(p_apex.squared_length());
  //NT r=1, d=2;//d=std::sqrt();
  
  /* Mink Sum 2D example 
  Polytope P;
  P.push_back(Hyperplane(Point(3,2),Direction(0,-1)));
  P.push_back(Hyperplane(Point(3,2),Direction(-1,0)));
  P.push_back(Hyperplane(Point(-2,-3),Direction(1,0)));
  P.push_back(Hyperplane(Point(-2,-3),Direction(0,1)));
  P.push_back(Hyperplane(Point(0,-3),Direction(-1,1)));
  //sandwitch
  NT r=1, d=std::sqrt(13.0);
  /**/
  
  // Random walks in K_i := the intersection of the ball i with P
  // the number of random points to be generated in each K_i
  int rnum = std::pow(e,-2) * 400 * n * std::log(n);
  int walk_len =  wl_c * std::pow(n,4);
  
  //rnum = e;
  //walk_len =  wl_c;
  
  int num_of_exp=1;
  for(int i=0; i<num_of_exp; ++i){
    std::cout<<n<<" "
		         <<rnum<<" "
		         <<walk_len<<" "<<std::flush;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    vars var(rnum,n,walk_len,err,0,0,0,0,rng,get_snd_rand,urdist,urdist1);
    //double v1 = volume1(P,var,var,r,d);
    tstop = (double)clock()/(double)CLOCKS_PER_SEC;
    //double v2 = volume2(P,n,rnum,walk_len,err,rng,get_snd_rand,urdist,urdist1);
  
	  //double exactvol = std::pow(2,n);
    double exactvol = std::pow(2,n)*std::pow(n,n)/factorial(n);

	  /*
	  std::cout<<rnum<<"\n\n\nALGORITHM 1\n-----------\nvolume = "
	           <<(1-e)*exactvol<<" < "<<v1<<" < "<<(1+e)*exactvol<<std::endl;
		std::cout<<"exact volume = "<<exactvol<<std::endl;
		std::cout<<"# walk steps = "<<walk_len<<std::endl;
		std::cout<<"# rand points = "<<rnum<<std::endl;
		std::cout<<"time = "<<tstop-tstart<<std::endl;
		*/
		std::cout.precision(15);
		std::cout
		         //<<v1<<" "
		         <<(1+e)*exactvol<<" "
		         <<exactvol<<" "
		         <<tstop-tstart<<std::endl;
		         
	}  
	
	
  //Minkowski_sum_naive(P1,P2,P);
  //for(V_polytope::iterator pit=P.begin(); pit!=P.end(); ++pit)
		//std::cout<<*pit<<std::endl;
	
	//V_polytope C = Vcross(n,-1,n);
  
	std::ofstream polymakefile;
	polymakefile.open("volume.polymake");
	//print_polymake_volfile(C,polymakefile);
  std::cout<<P[0]<<std::endl;
	print_polymake_volfile2(P,polymakefile);
	system ("polymake volume.polymake");
	std::cout<<std::endl;
  //}
  return 0;
}


