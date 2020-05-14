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


double HPolytopeCPP::compute_volume(char* method, int walk_len, double epsilon, uint seed){
   //Parameter setup
   int n = HP.dimension();
   int rnum;
   NT C;
   NT ratio;
   int N;
   int W;

   RNGType rng(seed);
   // boost::normal_distribution<> rdist(0,1);
   boost::random::uniform_real_distribution<>(urdist);
   boost::random::uniform_real_distribution<> urdist1(-1,1);

   Hpolytope HP_copy;
   HP_copy.init(HP.dimension(), HP.get_mat(), HP.get_vec());
   NT vol;

   if (strcmp(method,"sequence_of_balls")==0){
      rnum = std::pow(epsilon,-2) * 400 * n * std::log(n);
      C=2;
      ratio = 1.0-1.0/(NT(n));
      N = 500 * ((int) C) + ((int) (n * n / 2));
      W = 4*n*n+500;
      vars<NT, RNGType> var1_SoB(rnum,n,walk_len,1,0,1,0,0,0,0,rng,
         urdist,urdist1,-1.0,false,false,false,false,false,false,true,false);
      vol = volume(HP_copy, var1_SoB, CheBall);
   }
   else{//method=="gaussian_annealing"
      rnum = std::pow(1,-2) * 400 * n * std::log(n);
      C=2;
      ratio = 1.0-1.0/(NT(n));
      N = 500 * ((int) C) + ((int) (n * n / 2));
      W = 4*n*n+500;
      vars<NT, RNGType> var1_GA(rnum,n,walk_len,1,0,1,0,0,0,0,rng,
         urdist,urdist1,-1.0,false,false,false,false,false,false,true,false);
      vars_g<NT, RNGType> var2_GA(n,walk_len,N,W,1,epsilon,CheBall.second,rng,C,0.1,ratio,-1,false,
         false,false,false,false,false,false,true,false);
      vol = volume_cooling_gaussians(HP_copy, var2_GA, var1_GA, CheBall);
   }
   return (double)vol;
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
   double a_dummy = 1.0; //make this a parameter once gaussian_samples if a parameter too and can be true also.

   sampling_only<Point>(randPoints, HP, walk_len, n_samples, gaussian_samples, a_dummy, CheBall.first, var1, var2);

   auto n_si=0;
   for (auto it_s = randPoints.begin(); it_s != randPoints.end(); it_s++){
      for (auto i = 0; i != it_s->dimension(); i++){
         samples[n_si++] = (*it_s)[i];
      }
   }
}
