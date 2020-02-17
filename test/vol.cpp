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
#define VOLESTI_DEBUG
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "rotating.h"
#include "misc.h"
#include "linear_extensions.h"
#include "cooling_balls.h"
#include "cooling_hpoly.h"
#include "sample_only.h"
#include "exact_vols.h"

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

template <typename FT>
FT factorial(FT n)
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
    typedef double                    NT;
    typedef Cartesian<NT>             Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937            RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> Zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    int n, nexp=1, n_threads=1, W;
    int walk_len,N, nsam = 100, nu = 10, NNu;
    NT e=0.1;
    NT exactvol(-1.0), a=0.5;
    bool verbose=false, 
	 rand_only=false, 
	 round_only=false,
	 file=false, 
	 round=false, 
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
	 CB = false,
         birk=false,
         rotate=false,
         ball_walk=false,
         ball_rad=false,
         experiments=true,
         CG = false,
         Vpoly=false,
         Zono=false,
         cdhr=false,
         rdhr=false,
         user_randwalk = false,
         exact_zono = false,
         gaussian_sam = false,
         hpoly = false,
         billiard=false,
         win2 = false;

    //this is our polytope
    Hpolytope HP;
    Vpolytope VP; // RNGType only needed for the construction of the inner ball which needs randomization
    Zonotope  ZP;

    // parameters of CV algorithm
    bool user_W=false, user_N=false, user_ratio=false, user_NN = false, set_algo = false, set_error = false;
    NT ball_radius=0.0, diameter = -1.0, lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2, round_val = 1.0;
    NT C=2.0,ratio,frac=0.1,delta=-1.0,error=0.1;

  if(argc<2){
    std::cout<<"Use -h for help"<<std::endl;
    exit(-2);
  }
  
  //parse command line input vars
  for(int i=1;i<argc;++i){
      bool correct = false;

      if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")){
          std::cerr<<
                      "Usage:\n"<<
                      "-v, --verbose \n"<<
                      "-rot : does only rotating to the polytope\n"<<
                      "-ro, --round_only : does only rounding to the polytope\n"<<
                      "-rand, --rand_only : only samples from the convex polytope\n"<<
                      "-nsample : the number of points to sample in rand_olny mode\n"<<
                      "-gaussian : sample with spherical gaussian target distribution in rand_only mode\n"<<
                      "-variance : the variance of the spherical distribution in spherical mode (default 1)\n"<<
                      "-f1, --file1 [filename_type_Ax<=b] [epsilon] [walk_length] [threads] [num_of_experiments], for H-polytopes\n"<<
                      "-f2, --file2 [filename_type_V] [epsilon] [walk_length] [threads] [num_of_experiments], for V-polytopes\n"<<
                      "-f3, --file3 [filename_type_G] [epsilon] [walk_length] [threads] [num_of_experiments], for Zonotopes\n"<<
                      "-fle, --filele : counting linear extensions of a poset\n"<<
                      //"-c, --cube [dimension] [epsilon] [walk length] [threads] [num_of_experiments]\n"<<
                      "-exact_vol : the exact volume. Only fo zonotopes\n"<<
                      "-r, --round : enables rounding of the polytope as a preprocess\n"<<
                      "-e, --error epsilon : the goal error of approximation\n"<<
                      "-w, --walk_len [walk_len] : the random walk length. The default value is 1 for CB and CG and 10+d/10 for SOB\n"<<
                      "-exp [#exps] : number of experiments (default 1)\n"<<
                      "-cg : use Cooling Gaussians algorithm for volume computation. This is the default choice for H-polytopes in dimension >200\n"<<
                      "-cb : use Cooling Bodies algorithm for volume computation. This is the default choice for V-polytopes, Zonotopes and H-polytopes in dimension <=200\n"<<
                      "-sob : use Sequence of Balls algorithm for volume computation\n"<<
                      "-w, --walk_len [walk_len] : the random walk length (default 1)\n"<<
                      "-rdhr : use random directions HnR\n"<<
                      "-cdhr : use coordinate directions HnR. This is the default choice for H-polytopes\n"<<
                      "-biw : use Billiard walk. This is the default choice for V-polytopes and zonotopes\n"<<
                      "-L : the maximum length of the billiard trajectory\n"<<
                      "-baw : use ball walk\n"<<
                      "-bwr : the radius of the ball walk (default r*chebychev_radius/sqrt(max(1.0, a_i)*dimension\n"<<
                      "-e, --error epsilon : the goal error of approximation\n"<<
                      "-win_len : the size of the open window for the ratios convergence (for CB and CG algorithms)\n"<<
                      "-C : a constant for the upper boud of variance/mean^2 in schedule annealing. The default values is 2 (for CG algorithm)\n"
                      "-N : the number of points to sample in each step of schedule annealing. The default value N = 500*C + dimension^2/2 (for CG algorithm)\n"<<
                      "-frac : the fraction of the total error to spend in the first gaussian. The default frac=0.1 (for CG algo)\n"<<
                      "-ratio : parameter of schedule annealing, larger ratio means larger steps in schedule annealing. The default 1-1/dimension (for CG algorithm)\n"<<
                      "-lb : lower bound for the volume ratios in CB algorithm. The default values is 0.1\n"<<
                      "-ub : upper bound for the volume ratios in CB algorithm. The default values is 0.15\n"<<
                      "-alpha : the significance level for the statistical test in CB algorithm\n"<<
                      "-nu : the degrees of freedom of t-student to use in the t-tests in CB algorithm. The default value is 10\n"<<
                      "-nuN : the degrees of freedom of t-student to use in the t-tests and the number of samples to perform the statistical tests in CB algorithm. The default values is 10 and 125 (when -biw is used) or 120 + d*d/10 (when other random walks are used)\n"<<
                      std::endl;
          return 0;
      }
      if(!strcmp(argv[i],"--cube")){
          exactvol = std::pow(2,n);
          correct = true;
      }
      if(!strcmp(argv[i],"--exact")){
          exactvol = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-exact_vol")){
          exact_zono = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")){
          verbose = true;
          std::cout<<"Verbose mode\n";
          correct = true;
      }
      if(!strcmp(argv[i],"-rand")||!strcmp(argv[i],"--rand_only")){
          rand_only = true;
          std::cout<<"Generate random points only\n";
          correct = true;
      }
      if(!strcmp(argv[i],"-rdhr")){
          user_randwalk = true;
          rdhr = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-cdhr")){
          user_randwalk = true;
          cdhr = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-biw")){
          user_randwalk = true;
          billiard = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-baw")){
          user_randwalk = true;
          ball_walk = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-bwr")){
          delta = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-hpoly")){
          hpoly = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-win_len")){
          W = atof(argv[++i]);
          user_W = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-ratio")){
          ratio = atof(argv[++i]);
          user_ratio = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-frac")){
          frac = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-C")){
          C = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-N_an")){
          N = atof(argv[++i]);
          user_N = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-nuN")){
          NNu = atof(argv[++i]);
          nu = atof(argv[++i]);
          user_NN = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-nsample")){
          nsam = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-gaussian")){
          gaussian_sam = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-variance")){
          a = atof(argv[++i]);
          a = 1.0 / (2.0 * a);
          correct = true;
      }
      //reading from file
      if(!strcmp(argv[i],"-f1")||!strcmp(argv[i],"--file1")){
          file = true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          std::vector<std::vector<NT> > Pin;
          inp.open(argv[++i],std::ifstream::in);
          read_pointset(inp,Pin);
          n = Pin[0][1]-1;
          HP.init(Pin);
          if (verbose && HP.num_of_hyperplanes()<100){
              std::cout<<"Input polytope: "<<n<<std::endl;
              HP.print();
          }
          correct = true;
      }
      if(!strcmp(argv[i],"-f2")||!strcmp(argv[i],"--file2")){
          file = true;
          Vpoly = true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          std::vector<std::vector<NT> > Pin;
          inp.open(argv[++i],std::ifstream::in);
          read_pointset(inp,Pin);
          //std::cout<<"d="<<Pin[0][1]<<std::endl;
          n = Pin[0][1]-1;
          VP.init(Pin);
          if (verbose && VP.num_of_vertices()<100){
              std::cout<<"Input polytope: "<<n<<std::endl;
              VP.print();
          }
          correct = true;
      }
      if(!strcmp(argv[i],"-f3")||!strcmp(argv[i],"--file3")){
          file = true;
          Zono = true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          std::vector<std::vector<NT> > Pin;
          inp.open(argv[++i],std::ifstream::in);
          read_pointset(inp,Pin);
          //std::cout<<"d="<<Pin[0][1]<<std::endl;
          n = Pin[0][1]-1;
          ZP.init(Pin);
          correct = true;
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
          file = true;
          std::cout<<"Reading input from file..."<<std::endl;
          std::ifstream inp;
          inp.open(argv[++i],std::ifstream::in);
          std::ofstream os ("order_polytope.ine",std::ofstream::out);
          linear_extensions_to_order_polytope(inp,os);

          std::ifstream inp2;
          inp2.open("order_polytope.ine",std::ifstream::in);
          std::vector<std::vector<NT> > Pin;
          read_pointset(inp2,Pin);
          n = Pin[0][1]-1;
          HP.init(Pin);
          std::cout<<"Input polytope: "<<n<<std::endl;
          linear_extensions = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-r")||!strcmp(argv[i],"--round")){
          round = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-e")||!strcmp(argv[i],"--error")){
          e = atof(argv[++i]);
          error = e;
          set_error = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-w")||!strcmp(argv[i],"-walk_len")){
          walk_len = atof(argv[++i]);
          user_walk_len = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-exp")){
          nexp = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-L")){
          diameter = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-lb")){
          lb = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-ub")){
          ub = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-prob")){
          p = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-alpha")){
          alpha = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-nu")){
          nu = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-t")||!strcmp(argv[i],"--threads")){
          n_threads = atof(argv[++i]);
          correct = true;
      }
      if(!strcmp(argv[i],"-NN")){
          std::cout<<"flann software is needed for this option. Experimental feature."
                  <<"Currently under development."<<std::endl;
          correct = true;
      }
      if(!strcmp(argv[i],"-ro")){
          round_only = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-birk_sym")){
          birk = true;
          correct = true;
      }
      //rotate the polytope randomly
      if(!strcmp(argv[i],"-rot")){
          rotate = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-cg")){
          CG = true;
          set_algo = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-cb")){
          CB = true;
          set_algo = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-sob")){
          set_algo = true;
          correct = true;
      }
      if(correct==false){
          std::cerr<<"unknown parameters \'"<<argv[i]<<
                     "\', try "<<argv[0]<<" --help"<<std::endl;
          exit(-2);
      }
      
  }

  if (exact_zono) {
      if (!Zono) {
          std::cout<<"Exact volume computation is available only for zonotopes."<<std::endl;
          exit(-1);
      }
      NT vol_ex = exact_zonotope_vol<NT>(ZP);
      std::cout<<"Zonotope's exact volume = "<<vol_ex<<std::endl;
      return 0;
  }

  if (!set_algo) {
      if (Zono || Vpoly) {
          CB = true;
      } else {
          if (n <= 200) {
              CB = true;
          } else {
              CG = true;
          }
      }
  } else {
      if (!CB && !CG) {
          if (!set_error) {
              e = 1.0;
              error = 1.0;
          }
      }
  }

  if (!user_randwalk) {
      if (Zono || Vpoly) {
          if (CB) {
              billiard = true;
          } else {
              rdhr = true;
          }
      } else {
          cdhr = true;
      }
  } else if (!CB && !CG && billiard) {
      std::cout<<"Billiard is not supported for SOB algorithm. CDHR is used."<<std::endl;
      billiard = false;
      cdhr = true;
  } else if (CG && billiard) {
      billiard = false;
      if (Zono || Vpoly) {
          std::cout<<"Billiard is not supported for CG algorithm. RDHR is used."<<std::endl;
          rdhr = true;
      } else {
          std::cout<<"Billiard is not supported for CG algorithm. CDHR is used."<<std::endl;
          cdhr = true;
      }
  }

  /* RANDOM NUMBERS */
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  RNGType rng(seed);
  boost::normal_distribution<> rdist(0,1);
  boost::random::uniform_real_distribution<>(urdist);
  boost::random::uniform_real_distribution<> urdist1(-1,1);

  //Compute chebychev ball//
  std::pair<Point, NT> InnerBall;
  double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
  if (Zono) {
      InnerBall = ZP.ComputeInnerBall();
      if(billiard && diameter < 0.0){
          ZP.comp_diam(diameter, 0.0);
      }
  } else if(!Vpoly) {
      InnerBall = HP.ComputeInnerBall();
      if (InnerBall.second < 0.0) {
          std::cout<<"Polytope is unbounded!"<<std::endl;
          return -1.0;
      }
      if(billiard && diameter < 0.0){
          HP.comp_diam(diameter, InnerBall.second);
      }
  }else{
      if(CB) {
          if (round) {

              //std::cout<<"rounding is on"<<std::endl;
              InnerBall.first = VP.get_mean_of_vertices();
              InnerBall.second = 0.0;
              vars <NT, RNGType> var2(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second,
                                      2 * VP.get_max_vert_norm(), rng, urdist, urdist1, -1, verbose, rand_only, round,
                                      NN, birk, false, false, true,false);
              std::pair <NT, NT> res_round = rounding_min_ellipsoid(VP, InnerBall, var2);
              round_val = res_round.first;

              round = false;
              InnerBall.second = 0.0;
              InnerBall.first = Point(n);
              get_vpoly_center(VP);
              rmax = VP.get_max_vert_norm();
              if(billiard && diameter < 0.0) {
                  VP.comp_diam(diameter, 0.0);
              }

          } else {
              InnerBall.second = 0.0;
              InnerBall.first = Point(n);
              get_vpoly_center(VP);
              rmax = VP.get_max_vert_norm();
              if(billiard &&  diameter < 0.0){
                  VP.comp_diam(diameter, 0.0);
              }
          }
      } else {
          InnerBall = VP.ComputeInnerBall();
          if(billiard && diameter < 0.0){
              VP.comp_diam(diameter, 0.0);
          }
      }
  }

  if (ball_walk && delta < 0.0) {
      if (CG || gaussian_sam) {
          delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
      } else {
          delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
      }
  }

  double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
  if(verbose) std::cout << "Inner ball time: " << tstop1 - tstart1 << std::endl;
  if(verbose){
      std::cout<<"Inner ball center is: "<<std::endl;
      for(unsigned int i=0; i<n; i++){
          std::cout<<InnerBall.first[i]<<" ";
      }
      std::cout<<"\nradius is: "<<InnerBall.second<<std::endl;
  }
  
  // Set the number of random walk steps

  if(!user_walk_len) {
      if(!CG && !CB) {
          walk_len = 10 + n / 10;
      }else{
          walk_len = 1;
      }
  }
  if(!user_NN) {
      if(billiard) {
          NNu = 125;
      } else {
          NNu = 120 + (n * n) / 10;
      }
  }
  if(!user_N)
      N = 500 * ((int) C) + ((int) (n * n / 2));
  if(!user_ratio)
      ratio = 1.0-1.0/(NT(n));
  if(!user_W){
      if (CB) {
          if (billiard) {
              W = 170;
          } else {
              W = 3 * n * n + 400;
          }
      } else if (CG) {
          W = 4 * n * n + 500;
      }
  }

  // Timings
  double tstart, tstop;

  /* CONSTANTS */
  //error in hit-and-run bisection of P 
  const NT err=0.0000000001;
  const NT err_opt=0.01;

  //bounds for the cube	
  const int lw=0, up=10000, R=up-lw;

  // If no file specified construct a default polytope
  if(!file){
      std::cout << "A file has to be given as input!\n" << std::endl;
      return -1.0;
  }

  // If rotate flag is on rotate the polytope
  if(rotate) {
      if (Zono) {
          rotating<MT>(ZP);
          std::cout << "\n--------------\nRotated Zonotope\n" << std::endl;
          ZP.print();
      } else if (!Vpoly) {
          rotating<MT>(HP);
          std::cout << "\n--------------\nRotated polytope\nH-representation\nbegin\n" << std::endl;
          HP.print();
      } else {
          rotating<MT>(VP);
          std::cout << "\n--------------\nRotated polytope\nV-representation\nbegin\n" << std::endl;
          VP.print();
      }
      return 0;
  }
  if (rand_only) {
      std::list <Point> randPoints;

      vars<NT, RNGType> var1(0, n, walk_len, 1, 0, 0, 0, 0.0, 0, InnerBall.second, diameter, rng,
                urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr,billiard);
      vars_g<NT, RNGType> var2(n, walk_len, N, W, 1, 0, InnerBall.second, rng, C, frac, ratio, delta,
              verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);

      double tstart11 = (double)clock()/(double)CLOCKS_PER_SEC;
      if (Zono) {
          sampling_only<Point>(randPoints, ZP, walk_len, nsam, gaussian_sam, a, InnerBall.first, var1, var2);
      } else if (!Vpoly) {
          sampling_only<Point>(randPoints, HP, walk_len, nsam, gaussian_sam, a, InnerBall.first, var1, var2);
      } else {
          sampling_only<Point>(randPoints, VP, walk_len, nsam, gaussian_sam, a, InnerBall.first, var1, var2);
      }
      double tstop11 = (double)clock()/(double)CLOCKS_PER_SEC;
      if(verbose) std::cout << "Sampling time: " << tstop11 - tstart11 << std::endl;
      return 0;
  }

  // the number of random points to be generated in each K_i
  int rnum = std::pow(e,-2) * 400 * n * std::log(n);
  
  //RUN EXPERIMENTS
  int num_of_exp=nexp;
  double sum_time=0;
  NT min,max,sum=0;
  std::vector<NT> vs;
  NT average, std_dev;
  double Chebtime, sum_Chebtime=double(0);
  NT vol;

  for(unsigned int i=0; i<num_of_exp; ++i){
      std::cout<<"Experiment "<<i+1<<" ";
      tstart = (double)clock()/(double)CLOCKS_PER_SEC;

      // Setup the parameters
      vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0.0,0,InnerBall.second,diameter,rng,
               urdist,urdist1,delta,verbose,rand_only,round,NN,birk,ball_walk,cdhr,rdhr,billiard);

      if(round_only) {
          // Round the polytope and exit
          std::pair <NT, NT> res_round;
          if (Zono){
              res_round = rounding_min_ellipsoid(ZP, InnerBall, var);
              std::cout << "\n--------------\nRounded Zonotpe\n" << std::endl;
              ZP.print();
          } else if (!Vpoly) {
              res_round = rounding_min_ellipsoid(HP, InnerBall, var);
              std::cout << "\n--------------\nRounded polytope\nH-representation\nbegin\n" << std::endl;
              HP.print();
          } else {
              res_round = rounding_min_ellipsoid(VP, InnerBall, var);
              std::cout << "\n--------------\nRounded polytope\nV-representation\nbegin\n" << std::endl;
              VP.print();
          }
          std::cout << "end\n--------------\n" << std::endl;
          return 0;
      } else {
          // Estimate the volume
          if (CG) {

              // setup the parameters
              vars <NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, err, e, 0, 0.0, 0, InnerBall.second, diameter, rng,
                                      urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr,
                                      rdhr,billiard);

              vars_g <NT, RNGType> var1(n, walk_len, N, W, 1, error, InnerBall.second, rng, C, frac, ratio, delta,
                                        verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);

              if (Zono) {
                  vol = volume_gaussian_annealing(ZP, var1, var2, InnerBall);
              } else if (!Vpoly) {
                  vol = volume_gaussian_annealing(HP, var1, var2, InnerBall);
              } else {
                  vol = volume_gaussian_annealing(VP, var1, var2, InnerBall);
              }

          } else if (CB) {
              vars_ban <NT> var_ban(lb, ub, p, rmax, alpha, W, NNu, nu, win2);
              if (Zono) {
                  if (!hpoly) {
                      vol = vol_cooling_balls(ZP, var, var_ban, InnerBall);
                  } else {
                      vars_g <NT, RNGType> varg(n, 1, 500 * 2.0 +  NT(n * n) / 2.0, 6 * n * n + 500, 1, e, InnerBall.second, rng, C, frac, ratio, delta,
                                                verbose, rand_only, false, false, birk, false, true, false);
                      vol = vol_cooling_hpoly < HPolytope < Point > > (ZP, var, var_ban, varg, InnerBall);
                  }
              } else if (!Vpoly) {
                  vol = vol_cooling_balls(HP, var, var_ban, InnerBall);
              } else {
                  vol = vol_cooling_balls(VP, var, var_ban, InnerBall);
              }
              if (vol < 0.0) {
                  throw "Simulated annealing failed! Try to increase the walk length.";
                  return vol;
              }
          } else {
              if (Zono) {
                  vol = volume(ZP, var, InnerBall);
              } else if (!Vpoly) {
                  vol = volume(HP, var, InnerBall);
              } else {
                  vol = volume(VP, var, InnerBall);
              }
          }
      }

      NT v1 = vol;

      tstop = (double)clock()/(double)CLOCKS_PER_SEC;

      // Statistics
      sum+=v1;
      if(i==0){max=v1;min=v1;}
      if(v1>max) max=v1;
      if(v1<min) min=v1;
      vs.push_back(v1);
      sum_time +=  tstop-tstart;
      sum_Chebtime += Chebtime;

      if(round)
          std::cout<<" (rounding is ON)";
      std::cout<<std::endl;

      //Compute Statistics
      average=sum/(i+1);
      std_dev=0;
      for(std::vector<NT>::iterator vit=vs.begin(); vit!=vs.end(); ++vit){
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
                  <<HP.num_of_hyperplanes()<<" "
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
                 <<HP.num_of_hyperplanes()<<" "
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
