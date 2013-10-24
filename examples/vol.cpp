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

#include <vol_rand.h>

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
	int n,nexp;
	double wl_c, e;
	//parse command line input
	bool verbose=false;
	bool file=false;
	stdHPolytope<double> P;
	
	for(int i=1;i<argc;++i){
		bool correct=false;
    if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")){
      std::cerr<<
        "Usage:\n"<<
        "-i [dimension] [epsilon] [walk length] [num of experiments]\n"<<
        "-v, --verbose \n"<<
        "-f, --file [filename] [epsilon] [walk length] [num of experiments]\n";
      exit(-1);
    }
		if(!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")){
      verbose=true;
      ++i;
      std::cout<<"Verbose mode\n";
      correct=true;
    } 
    if(!strcmp(argv[i],"-f")||!strcmp(argv[i],"--file")){
      if(argc-i<4){
	      std::cerr<<"Wrong number of arguments \'"<<argc-i-1<<
	        "\', should be 3,"<<" Try also --help"<<std::endl;
	      exit(-2);
			}else{
	      file=true;
	      std::ifstream inp;
	      std::vector<std::vector<double> > Pin;
	      inp.open(argv[++i],std::ifstream::in);
	      read_pointset(inp,Pin);
	      //std::cout<<"d="<<Pin[0][1]<<std::endl;
	      n = Pin[0][1]-1;
	      P.init(Pin);
	      P.print();
	      //constants
				e = atof(argv[++i]);
				wl_c = atof(argv[++i]);
				nexp = atof(argv[++i]);
	      //exit(1);
	      correct=true;
	    }
    }
    if(!strcmp(argv[i],"-i")||!strcmp(argv[i],"--input")){
      if(argc-i<5){
				std::cerr<<"Wrong number of arguments \'"<<argc-i-1<<
        "\', should be 4,"<<" Try also --help"<<std::endl;
      exit(-2);
			}else{
	      //dimension
				n = atoi(argv[++i]);
				//constants
				e = atof(argv[++i]);
				wl_c = atof(argv[++i]);
				nexp = atof(argv[++i]);
				correct=true;
			}  
    }
    if(correct==false){
      std::cerr<<"unknown parameter \'"<<argv[i]<<
        "\', try "<<argv[0]<<" --help"<<std::endl;
      exit(-2);
    }
		
	}
	

  //test zone
  /*
  stdHPolytope<double> P2(10);
  P2.print();
  std::vector<double> p2;
  //for (int i=0; i<10; ++i)
  p2.push_back(double(1));
  p2.push_back(double(-1));
  p2.push_back(double(0.99));
  p2.push_back(double(1));
  p2.push_back(double(-1));
  p2.push_back(double(-1));
  p2.push_back(double(1));
  p2.push_back(double(-1));
  p2.push_back(double(-1));
  p2.push_back(double(-1));
  Point pp2(10,p2.begin(),p2.end());
  std::cout<<P2.is_in(pp2)<<std::endl;
  exit(1);
  /**/
	//timings
  double tstart, tstop;


  /* CONSTANTS */
  //error in hit-and-run bisection of P 
  const double err=0.0000000001; 
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

  //FOR EXPERIMENTS
	//for(n=2; n<13; n+=2){

	/* VOLUME */
  
  /* CUBE */
	//Polytope P = cube(n,-1,1);
	
  if(!file)
    P.init(n);
  //sandwitch
  std::vector<NT> coords_apex(n,1);
	Vector p_apex(n,coords_apex.begin(),coords_apex.end());
  NT r=1, d=std::sqrt(p_apex.squared_length());
  /**/
  
  /*  Cross-polytope 
  Polytope P=cross(n,-1,n);
  std::vector<NT> apex;
	for(int j=0; j<n; ++j){
		if(j==0) apex.push_back(NT(n));
		else apex.push_back(NT(0));
  }
  Vector p_apex(n,apex.begin(),apex.end());
  NT r=1, d=std::sqrt(p_apex.squared_length());
  //NT r=1, d=2;//d=std::sqrt();
  /**/
  
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
  //int rnum = e;
  
  // The number of hit-&-run steps applied to each point   
  //int walk_len =  wl_c * std::pow(n,4);
  int walk_len =  wl_c;
  
  //RUN EXPERIMENTS
  int num_of_exp=nexp;
  double sum=0, sum_time=0;
  double min,max;
  std::vector<double> vs;
  double average, std_dev, exactvol;
  
  for(int i=0; i<num_of_exp; ++i){
    std::cout<<"Experiment "<<i+1<<" ";
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    vars var(rnum,n,walk_len,err,0,0,0,0,rng,get_snd_rand,urdist,urdist1,verbose);
    double v1 = volume1_reuse(P,var,var,r,d);
    tstop = (double)clock()/(double)CLOCKS_PER_SEC;
    //double v2 = volume2(P,n,rnum,walk_len,err,rng,get_snd_rand,urdist,urdist1);
     
	  //Used to Compute Statistics
    sum+=v1;
    if(i==0){max=v1;min=v1;}
    if(v1>max) max=v1;
    if(v1<min) min=v1;
    vs.push_back(v1);
		sum_time +=  tstop-tstart;
		std::cout<<"\t vol= "<<v1<<"\t time= "<<tstop-tstart<<std::endl;        
	  
		//Compute Statistics
		average=sum/(i+1);
		std_dev=0;
		for(std::vector<double>::iterator vit=vs.begin(); vit!=vs.end(); ++vit){
			std_dev += std::pow(*vit - average,2);
		}
		std_dev = std::sqrt(std_dev/(i+1));
		
		exactvol = std::pow(2,n);
	  //double exactvol = std::pow(2,n)*std::pow(n,n)/factorial(n);
		std::cout.precision(7);
		
		//Print statistics
		std::cout<<"STATISTICS:"<<std::endl;
		//std::cout<<"d #experiments (1-e)vol, (1+e)vol vol N walklen average min, max std_dev (vol-v*)/vol t"<<std::endl;
		std::cout 
		           <<n<<" "
		           <<num_of_exp<<" "
		           <<exactvol<<" ["
		           <<(1-e)*exactvol<<", "
			         <<(1+e)*exactvol<<"] "
		           <<rnum<<" "
		           <<walk_len<<" "
			         <<average<<" ["
			         <<min<<", "
			         <<max<<"] "
			         <<std_dev<<" "
			         <<(exactvol-average)/exactvol<<" "
			         <<(max-min)/average<<" "
			         <<sum_time/(i+1)<<std::endl;
	}
	/*
  // EXACT COMPUTATION WITH POLYMAKE
  /*
	std::ofstream polymakefile;
	polymakefile.open("volume.polymake");
	//print_polymake_volfile(C,polymakefile);
  std::cout<<P[0]<<std::endl;
	print_polymake_volfile2(P,polymakefile);
	system ("polymake volume.polymake");
	std::cout<<std::endl;
  */
  //}
  
  return 0;
}


