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
  double tstart, tstop, tstart2, tstop2;

  /* CONSTANTS */
  // L = log(R/r)
  const int L=30;
  //error in hit-and-run bisection of P 
  const double err=0.001; 
  const double err_opt=0.00001; 
  const double err_opt_bisect=0.0049;
  //bounds for the cube	
  const int lw=0;//, R=up-lw;
  
		
		
  /* RANDOM NUMBERS */
  // the random engine with time as a seed
  RNGType rng((double)time(NULL));
  // standard normal distribution with mean of 0 and standard deviation of 1 
	boost::normal_distribution<> rdist(0,1); 
	boost::variate_generator< RNGType, boost::normal_distribution<> > 
											get_snd_rand(rng, rdist);
											
	//exponential distribution
	//boost::exponential_distribution<> rdist(0.1); //i.e. lambda
	//boost::variate_generator< RNGType, boost::exponential_distribution<> > 
	//										get_snd_rand(rng, rdist);	
																	 
  // uniform distribution 
  boost::random::uniform_real_distribution<>(urdist); 
  boost::random::uniform_real_distribution<> urdist1(-1,1); 
  
  //cross(n,-1,1);
  //exit(1);
  
  for (n=2; n<12; ++n){
  //n=11;		
		
  /* OPTIMIZATION */
  /*!!! given a direction w compute a vertex v of K that maximize w*v */
  
  //Define the input polytope and the optimization direction
  
  /*  Cube   
  int up=NT(100);
  Polytope P=cube(n,-1,1);
  std::vector<NT> ww(n,1);
  //w=w/w.squared_length();           //normalize
  Vector w(n,ww.begin(),ww.end());
  Vector opt = w;
  /**/
  /*  Cross-polytope */
  int up=NT(pow(2,n));
  Polytope P=cross_skinny2(n,-1,1);
  std::vector<NT> ww;
	for(int j=0; j<n; ++j){
		if(j==0) ww.push_back(NT(1));
		else ww.push_back(NT(0));
  }
  
  Vector w(n,ww.begin(),ww.end());
  
  std::vector<NT> optt;
	for(int j=0; j<n; ++j){
		if(j==0) optt.push_back(NT(pow(2,n)));
		else optt.push_back(NT(0));
  }
  Vector opt(n,optt.begin(),optt.end());
  /**/
  
  /* Mink Sum 2D polar example 
  int up=NT(1);
  Polytope P;
  P.push_back(Hyperplane(Point(0.0,0.5),Direction(1.0,-1.0)));
  P.push_back(Hyperplane(Point(0.0,0.5),Direction(-3.0,-2.0)));
  P.push_back(Hyperplane(Point(NT(1.0/3.0),NT(-1.0/3.0)),Direction(-1.0,0.0)));
  P.push_back(Hyperplane(Point(NT(1.0/3.0),NT(-1.0/3.0)),Direction(0.0,1.0)));
  P.push_back(Hyperplane(Point(-0.5,0.0),Direction(2.0,3.0)));
  Vector w(0.0,1.0);
  Vector opt(0.5,0.0);
  /**/
  
  //the feasible point that approximates opt (at the end fp should be opt)	
  Point fp1, fp2, fp3;
    
  //the number of generated random points 
  int rnum =  e * n * std::pow(std::log(n),2);;
  //the number of steps in every random walk
  
  int walk_len =  wl_c * std::pow(n,4);
  double t1=0, t2=0, t3=0, max_err1=0, max_err2=0, max_err3=0;
  int num_of_exp = 1;
    
  vars var1(rnum,n,walk_len,err,err_opt_bisect,lw,up,L,rng,get_snd_rand,urdist,urdist1);
  vars var2(rnum,n,walk_len,err,err_opt,lw,up,L,rng,get_snd_rand,urdist,urdist1);
  vars var3(rnum,n,walk_len,err,err_opt,lw,up,L,rng,get_snd_rand,urdist,urdist1);

  std::vector<NT> test(n,0);
  //std::cout<<"feasible="<<Sep_Oracle(P,Point(n,test.begin(),test.end()),var1).get_is_in()<<std::endl;
 
  //for(Polytope::iterator it=P.begin(); it!=P.end(); ++it){
  //  std::cout<<(*it).has_on_positive_side(Point(n,test.begin(),test.end()))<<std::endl;
  //  std::cout<<(*it).orthogonal_direction()<<std::endl;
  //}
  // exit(1);
  
  //std::cout<<"# walk steps = "<<walk_len<<std::endl;
	//std::cout<<"# rand points = "<<rnum<<std::endl<<std::endl;
  //std::cout.precision(3);
  std::cout<<n<<"\t "
	           <<e<<"\t "
	           <<wl_c<<"\t "
		         <<rnum<<"\t "
		         <<walk_len<<"\t "<<std::flush;
  for(int i=0; i<num_of_exp; ++i){
	  //do optimization 
	  tstart = (double)clock()/(double)CLOCKS_PER_SEC;
	  //opt_bisect(P,var1,fp1,w);
	  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
	  t1 += tstop-tstart;
	  if (max_err1 < std::abs((fp1-CGAL::Origin())*w - opt*w))
	    max_err1=std::abs((fp1-CGAL::Origin())*w - opt*w);
	    
	  //std::cout<<"Interior point"<<std::endl;
	  //do interior point optimization	
		tstart = (double)clock()/(double)CLOCKS_PER_SEC;
	  optimization(P,var2,fp2,w);
	  //feasibility(P,fp2,var2);
	  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
	  t2 += tstop-tstart;
	  if (max_err2 < std::abs((fp2-CGAL::Origin())*w - opt*w))
	    max_err2=std::abs((fp2-CGAL::Origin())*w - opt*w);
	  //std::cout<<"t1="<<tstop-tstart<<" "<<t1<<" t2="<<tstop2-tstart2<<" "<<t2<<std::endl;
	  
	  tstart = (double)clock()/(double)CLOCKS_PER_SEC;
	  //opt_interior(P,var3,fp3,w);
	  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
	  t3 += tstop-tstart;
	  if (max_err3 < std::abs((fp3-CGAL::Origin())*w - opt*w))
	    max_err3=std::abs((fp3-CGAL::Origin())*w - opt*w);
	  //std::cout<<"t1="<<tstop-tstart<<" "<<t1<<" t2="<<tstop2-tstart2<<" "<<t2<<std::endl;
	  
  }
  /*
  //print the results
  std::cout<<"------------------"<<std::endl;

  std::cout<<"OPT = "<<fp<<std::endl;
  std::cout<<"max err = "<<max_err1<<std::endl;
  std::cout<<"time = "<<t1/num_of_exp<<std::endl<<std::endl;
	//std::cout<<"max num of iterations = "<<2*n*L<<std::endl;
  //
  std::cout<<"OPT I="<<fp2<<std::endl;
  std::cout<<"max err = "<<max_err2<<std::endl;
	std::cout<<"time = "<<t2/num_of_exp<<std::endl;
	*/
	std::cout  <<max_err1<<"\t "
		         <<t1/num_of_exp<<"\t "
		         <<max_err2<<"\t "
		         <<t2/num_of_exp<<"\t "
             <<max_err3<<"\t "
		         <<t3/num_of_exp<<std::endl;
  }
  /**/
  
  //std::vector<NT> testp(n,NT(0.2));
  //std::cout<<B.has_on_positive_side(Point(n,testp.begin(),testp.end()))<<std::endl;
  
  
	  
  return 0;
}


