#include <iostream>
#include <math.h>
#include "bindings.h"

using namespace std;


HPolytopeCPP::HPolytopeCPP() {}
HPolytopeCPP::HPolytopeCPP(double *A_np, double *b_np, int n_hyperplanes, int n_variables){
   MT A;
   VT b;
   A.resize(n_hyperplanes,n_variables);
   b.resize(n_hyperplanes);

   int index=0;
   for (int i=0; i<n_hyperplanes; i++){
      b(i) = b_np[i];
      for (int j=0; j<n_variables; j++){
         A(i,j) = A_np[index];
         index++;
      }
   }

   HP.init(n_variables,A,b);
   CheBall = HP.ComputeInnerBall();
}

HPolytopeCPP::~HPolytopeCPP(){}


double HPolytopeCPP::compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, uint seed){
  
// the following command used to be like this "<boost::mt19937, NT, 3>" but we removed "3"   
   typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
   double volume;
// strcmp returns a lexical difference (short-circuit serial byte comparator) of the two strings you have given as parameters. 0 means that both strings are equal  
// we are about to have 3 methods for computing the volume () and a number of random walks for each of those
   if (strcmp(vol_method,"sequence_of_balls")==0)){   // we have 
      
      if (strcmp(walk_method,"uniform_ball"==0)){
         volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"CDHR"==0)){
         volume = volume_sequence_of_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"RDHR"==0)){
         volume = volume_sequence_of_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
         
      }
   }
   else if (strcmp(vol_method,"cooling_gaussian")==0){
      if (strcmp(walk_method,"gaussian_ball"==0){
         volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_CDHR"==0){
         volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_RDHR"==0){
         volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   
   } else if (strcmp(vol_method,"cooling_balls")==0){
       
       if (strcmp(walk_method,"uniform_ball"==0){
         volume = volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len);
       } else if (strcmp(walk_method,"CDHR"==0){
         volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, e, walk_len);
       } else if (strcmp(walk_method,"RDHR"==0){
         volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, e, walk_len);
       } else if (strcmp(walk_method,"billiard"==0){
         volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len);
       }
       
      
      
      
      
      
   }
   return volume;
}












void HPolytopeCPP::generate_samples(int walk_len, int n_samples, double* samples, uint seed){
   //Parameter setup
   int n = HP.dimension();

   int rnum = std::pow(1,-2) * 400 * n * std::log(n);
   NT C=2;
   NT ratio = 1.0-1.0/(NT(n));
   int N = 500 * ((int) C) + ((int) (n * n / 2));
   int W = 4*n*n+500;

   double epsilon_dummy = 0.1;

   RNGType rng(seed);
   // boost::normal_distribution<> rdist(0,1);
   boost::random::uniform_real_distribution<>(urdist);
   boost::random::uniform_real_distribution<> urdist1(-1,1);


   vars<NT, RNGType> var1(rnum,n,walk_len,1,0,1,0,0,0,CheBall.second,rng,
      urdist,urdist1,-1,false,false,false,false,false,false,true,false);
   vars_g<NT, RNGType> var2(n,walk_len,N,W,1,epsilon_dummy,CheBall.second,rng,
      C,0.1,ratio,-1,false,false,false,false,false,false,false,true,false);

   std::list <Point> randPoints;
   bool gaussian_samples = false;
   //make this a parameter once gaussian_samples if a parameter too and can be true also.
   double a_dummy = 1.0;

   uniform_sampling<Point>(randPoints, HP, walk_len, n_samples, gaussian_samples,
                           a_dummy, CheBall.first, var1, var2);


// The following block of code should NOT be removed!
   auto n_si=0;
   for (auto it_s = randPoints.begin(); it_s != randPoints.end(); it_s++){
      for (auto i = 0; i != it_s->dimension(); i++){
         samples[n_si++] = (*it_s)[i];
      }
   }
}
