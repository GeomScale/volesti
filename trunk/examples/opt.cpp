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
	/* VARIABLES */
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
	
	//for timing
  double tstart, tstop;

  /* CONSTANTS */
  // L = log(R/r)
  const int L=30;
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
  
  
  
  /* OPTIMIZATION */
  /*!!! given a direction w compute a vertex v of K that maximize w*v */
  
  //this is the input polytope
  Polytope P=cube(n,0,10);
  
  //this is the optimization direction
  std::vector<NT> ww(n,1);
  Vector w(n,ww.begin(),ww.end());
  w=w/w.squared_length();           //normalize
  
  //the feasible point that approximates opt (at the end fp should be opt)
  Point fp;
  //the number of generated random points 
  int rnum =  e * n * std::pow(std::log(n),2);;
  //the number of steps in every random walk
  int walk_len =  wl_c * std::pow(n,4);
  
  //compute the optimum 
  tstart = (double)clock()/(double)CLOCKS_PER_SEC;
  optimization(P,rnum,n,walk_len,err,lw,up,L,rng,get_snd_rand,urdist,urdist1,fp,w);
  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
  
  //print the results
  std::cout<<"OPT="<<fp<<std::endl;
  std::cout<<"# walk steps = "<<walk_len<<std::endl;
	std::cout<<"# rand points = "<<rnum<<std::endl;
	std::cout<<"time = "<<tstop-tstart<<std::endl;
  /**/
  
  //std::vector<NT> testp(n,NT(0.2));
  //std::cout<<B.has_on_positive_side(Point(n,testp.begin(),testp.end()))<<std::endl;
  
  /* Optimization with bisection
	 * TODO: make it a function!!!
	 */
  
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


