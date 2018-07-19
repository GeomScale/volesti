// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#include "Eigen/Eigen"
#include "volume.h"

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
	//Deafault values
    int n, nexp=1, n_threads=1;
	int walk_len;//to be defined after n
    double e=1;
    	double exactvol(-1.0);
    bool verbose=false, 
	 rand_only=false, 
	 round_only=false,
	 file=false, 
	 round=false, 
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         experiments=true,
         coordinate=true;
	
	//this is our polytope
	Polytope<NT> P;
	int magnitude=0;
	bool exper=false;
	
  if(argc<2){
    std::cout<<"Use -h for help"<<std::endl;
    exit(-2);
  }
  
  //parse command line input vars
  for(int i=1;i<argc;++i){
      bool correct=false;

      if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")){
          std::cerr<<
                      "Usage:\n"<<
                      "-v, --verbose \n"<<
                      "-rdhr : use random directions HnR, default is coordinate directions HnR\n"
                      "-rand, --rand_only : generates only random points\n"<<
                      "-f1, --file1 [filename_type_Ax<=b] [epsilon] [walk_length] [threads] [num_of_experiments]\n"<<
                      //"-f2, --file2 [filename_type_Ax=b,x>=0] [epsilon] [walk_length] [threads] [num_of_experiments]\n"<<
                      "-fle, --filele : counting linear extensions of a poset\n"<<
                      //"-c, --cube [dimension] [epsilon] [walk length] [threads] [num_of_experiments]\n"<<
                      "--exact : the exact volume\n"<<
                      "--cube : input polytope is a cube\n"<<
                      "-r, --round : enables rounding of the polytope as a preprocess\n"<<
                      "-ro, --round_only : does only rounding to the polytope\n"<<
                      "-e, --error epsilon : the goal error of approximation\n"<<
                      "-w, --walk_len [walk_len] : the random walk length (default 10)\n"<<
                      "-exp [#exps] : number of experiments (default 1)\n"<<
                      "-t, --threads #threads : the number of threads to be used\n"<<
                      "-ΝΝ : use Nearest Neighbor search to compute the boundary oracles\n"<<
                      "-birk_sym : use symmetry to compute more random points (only for Birkhoff polytopes)\n"<<
                      std::endl;
          return 0;
      }
      if(!strcmp(argv[i],"--cube")){
          exactvol = std::pow(2,n);
          //exactvol = std::pow(2,n)/std::tgamma(n+1);//factorial of a natural number n is gamma(n+1)
          correct=true;
      }
      if(!strcmp(argv[i],"--exact")){
          exactvol = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-mag")){
          magnitude=int(atof(argv[++i]));
          exper=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")){
          verbose=true;
          std::cout<<"Verbose mode\n";
          correct=true;
      }
      if(!strcmp(argv[i],"-rand")||!strcmp(argv[i],"--rand_only")){
          rand_only=true;
          std::cout<<"Generate random points only\n";
          correct=true;
      }
      if(!strcmp(argv[i],"-rdhr")){
          coordinate=false;
          correct=true;
      }
      //reading from file
      if(!strcmp(argv[i],"-f1")||!strcmp(argv[i],"--file1")){
          file=true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          std::vector<std::vector<double> > Pin;
          inp.open(argv[++i],std::ifstream::in);
          read_pointset(inp,Pin);
          //std::cout<<"d="<<Pin[0][1]<<std::endl;
          n = Pin[0][1]-1;
          P.init(Pin);
          if (verbose && P.num_of_hyperplanes()<100){
              std::cout<<"Input polytope: "<<n<<std::endl;
              P.print();
          }
          correct=true;
      }
      /*
    if(!strcmp(argv[i],"-f2")||!strcmp(argv[i],"--file2")){
            file=true;
      std::ifstream inp;
      std::vector<std::vector<double> > Pin;
      inp.open(argv[++i],std::ifstream::in);
      read_pointset(inp,Pin);
      //std::cout<<"d="<<Pin[0][1]<<std::endl;
      //n = Pin[0][1]-1;
      P.init(Pin);
      P.rref();
      n=P.dimension();
      //if (verbose && P.num_of_hyperplanes()<1000){
    // std::cout<<"Input polytope: "<<n<<std::endl;
      //  P.print();
      //}
      correct=true;
    }
*/
      //reading linear extensions and order polytopes
      if(!strcmp(argv[i],"-fle")||!strcmp(argv[i],"--filele")){
          file=true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          inp.open(argv[++i],std::ifstream::in);
          std::ofstream os ("order_polytope.ine",std::ofstream::out);
          linear_extensions_to_order_polytope(inp,os);

          std::ifstream inp2;
          inp2.open("order_polytope.ine",std::ifstream::in);
          std::vector<std::vector<double> > Pin;
          read_pointset(inp2,Pin);

          //std::cout<<"d="<<Pin[0][1]<<std::endl;
          n = Pin[0][1]-1;
          P.init(Pin);
          //if (verbose && P.num_of_hyperplanes()<100){
          std::cout<<"Input polytope: "<<n<<std::endl;
          //P.print();
          //}
          linear_extensions = true;
          correct=true;
      }
      if(!strcmp(argv[i],"-r")||!strcmp(argv[i],"--round")){
          round = true;
          correct=true;
      }
      if(!strcmp(argv[i],"-e")||!strcmp(argv[i],"--error")){
          e = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-w")||!strcmp(argv[i],"--walk_len")){
          walk_len = atof(argv[++i]);
          user_walk_len=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-exp")){
          nexp = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-t")||!strcmp(argv[i],"--threads")){
          n_threads = atof(argv[++i]);
          correct=true;
      }
      if(!strcmp(argv[i],"-NN")){
          /*
      if (verbose) std::cout<<"Building search data-srtuctures..."<<std::endl;
      NN=true;
      P.dual(1);
      P.dual(-1);
      */
          std::cout<<"flann software is needed for this option. Experimental feature."
                  <<"Currently under development."<<std::endl;
          correct=true;
      }
      if(!strcmp(argv[i],"-ro")){
          round_only=true;
          correct=true;
      }
      if(!strcmp(argv[i],"-birk_sym")){
          birk=true;
          correct=true;
      }
      //rotate the polytope randomly
      if(!strcmp(argv[i],"-rot")){
          rotate=true;
          correct=true;
      }
      if(correct==false){
          std::cerr<<"unknown parameters \'"<<argv[i]<<
                     "\', try "<<argv[0]<<" --help"<<std::endl;
          exit(-2);
      }
      
  }//for i
  
  //Compute chebychev ball//
  double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
  std::pair<Point,double> CheBall = solveLP(P);
  double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
  if(verbose) std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;
  if(verbose){
      std::cout<<"Chebychev center is: "<<std::endl;
      for(int i=0; i<P.dimension(); i++){
          std::cout<<CheBall.first[i]<<" ";
      }
      std::cout<<"\nradius is: "<<CheBall.second<<std::endl;
  }
  
  // Set the number of random walk steps
  if(!user_walk_len)
      walk_len=10 + n/10;

  // Timings
  double tstart, tstop;

  /* CONSTANTS */
  //error in hit-and-run bisection of P 
  const double err=0.0000000001; 
  const double err_opt=0.01; 
  //bounds for the cube	
  const int lw=0, up=10000, R=up-lw;
  
   /* RANDOM NUMBERS */
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  //RNGType rng(std::time(0));
  RNGType rng(seed);
  boost::normal_distribution<> rdist(0,1);
  boost::random::uniform_real_distribution<>(urdist);
  boost::random::uniform_real_distribution<> urdist1(-1,1);

  // If no file specified construct a default polytope
  if(!file){
      P.init(n);
  }

  // If rotate flag is on rotate the polytope
  if(rotate){
      rotating(P);
      //P.print();
  }

  // Random walks in K_i := the intersection of the ball i with P
  // the number of random points to be generated in each K_i
  int rnum = std::pow(e,-2) * 400 * n * std::log(n);
  
  //RUN EXPERIMENTS
  int num_of_exp=nexp;
  double sum=0, sum_time=0;
  double min,max;
  std::vector<double> vs;
  double average, std_dev;
  double Chebtime, sum_Chebtime=double(0);
  NT vol;
  
  for(int i=0; i<num_of_exp; ++i){
      std::cout<<"Experiment "<<i+1<<" ";
      Polytope<double> P_to_test(P);
      tstart = (double)clock()/(double)CLOCKS_PER_SEC;

      // Setup the parameters
      vars var(rnum,n,walk_len,n_threads,err,e,0,0,0,rng,
               urdist,urdist1,verbose,rand_only,round,NN,birk,coordinate);

      if(round_only){
          // Round the polytope and exit
          std::pair<double,double> res_round;
          res_round = rounding_min_ellipsoid(P,CheBall.first,CheBall.second,var);
          double round_value = res_round.first;
          std::cout<<"\n--------------\nRounded polytope\nH-representation\nbegin\n"<<std::endl;
          P.print();
          std::cout<<"end\n--------------\n"<<std::endl;
      }else{
          // Estimate the volume
          vol = volume(P_to_test,var,var,CheBall);
          //if(rotate) vol = std::sqrt(vol);
          //std::cout<<vol<<std::endl;
      }

      double v1 = vol;

      tstop = (double)clock()/(double)CLOCKS_PER_SEC;
      //double v2 = volume2(P,n,rnum,walk_len,err,rng,get_snd_rand,urdist,urdist1);

      // Statistics
      sum+=v1;
      if(i==0){max=v1;min=v1;}
      if(v1>max) max=v1;
      if(v1<min) min=v1;
      vs.push_back(v1);
      sum_time +=  tstop-tstart;
      sum_Chebtime += Chebtime;

      //std::cout<<"\t vol= "<<v1<<"\t time= "<<tstop-tstart;
      if(round)
          std::cout<<" (rounding is ON)";
      std::cout<<std::endl;

      //Compute Statistics
      average=sum/(i+1);
      std_dev=0;
      for(std::vector<double>::iterator vit=vs.begin(); vit!=vs.end(); ++vit){
          std_dev += std::pow(*vit - average,2);
      }
      std_dev = std::sqrt(std_dev/(i+1));

      std::cout.precision(7);

      //MEMORY USAGE
      //struct proc_t usage;
      //look_up_our_self(&usage);

      //Print statistics
      //std::cout<<"\nSTATISTICS:"<<std::endl;
      if (!experiments){
          std::cout
                  <<"Dimension= "
                  <<n<<" "
                   //<<argv[]<<" "
                  <<"\nNumber of hyperplanes= "
                  <<P.num_of_hyperplanes()<<" "
                  <<"\nNumber of runs= "
                  <<num_of_exp<<" "
                  <<"\nError parameter= "
                  <<e
                  <<"\nTheoretical range of values= "<<" ["
                  <<(1-e)*exactvol<<","
                  <<(1+e)*exactvol<<"] "
                  <<"\nNumber of random points generated in each iteration= "
                  <<rnum<<" "
                  <<"\nRandom walk length= "
                  <<walk_len<<" "
                  <<"\nAverage volume (avg)= "
                  <<average
                  <<"\nmin,max= "
                    " ["
                  <<min<<","
                  <<max<<"] "
                  <<"\nStandard deviation= "
                  <<std_dev<<" "
                  <<"\n(max-min)/avg= "
                  <<(max-min)/average<<" "
                  <<"\nTime(sec)= "
                  <<sum_time/(i+1)<<" "
                  <<"\nTime(sec) Chebyshev= "
                  <<sum_Chebtime/(i+1)<<" "
                    //<<usage.vsize
                  <<std::endl;
    
      if(exactvol!=-1.0){
	      std::cout 
	           <<"\nExact volume= "
	           <<exactvol<<" "
	           <<"\n(vol-avg)/vol= "
	           <<(exactvol-average)/exactvol<<" "
               <<std::endl;
      }
	} else 
    	std::cout 
                 <<n<<" "
                 //<<argv[]<<" "
                 <<P.num_of_hyperplanes()<<" "
                 <<num_of_exp<<" "
                 <<exactvol<<" "
                 <<e<<" ["
                 <<(1-e)*exactvol<<","
                 <<(1+e)*exactvol<<"] "
                 <<rnum<<" "
                 <<walk_len<<" "
                 <<average<<" ["
                 <<min<<","
                 <<max<<"] "
                 <<std_dev<<" "
                 <<(exactvol-average)/exactvol<<" "
                 <<(max-min)/average<<" "
                 <<sum_time/(i+1)<<" "
                 <<sum_Chebtime/(i+1)<<" "
                 //<<usage.vsize
                 <<std::endl;
	}
	if(exper){
		
		if (int( std::ceil( -std::log10( average ) ) )<0){
			if( magnitude==0){
				std::cout<<"TEST PASSED"<<std::endl;
			}else{
				std::cout<<"TEST FAILED"<<std::endl;
			}
		}else if( magnitude==int( std::ceil( -std::log10( average ) ) ) ){
			std::cout<<"TEST PASSED"<<std::endl;
		}else{
			std::cout<<"TEST FAILED"<<std::endl;
		}
		std::cout<<"-----------------------------------------------\n";
		std::cout<<" # # # # # # # # # # # # # # # # # # # # # # #\n";
		std::cout<<"-----------------------------------------------\n";
	}
	
  if(linear_extensions)
		   std::cout <<"Number of linear extensions= "<<vol*factorial(n)<<std::endl;
  
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
