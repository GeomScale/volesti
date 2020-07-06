#include <iostream>
#include <math.h>
#include "bindings.h"

using namespace std;
//unsigned int dim;

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
   //dim = HP.dimension();
}
HPolytopeCPP::~HPolytopeCPP(){}

//////////          This is where the "COMPUTE_VOLUME" class starts  

double HPolytopeCPP::compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed){
  
   //typedef BoostRandomNumberGenerator<boost::mt19937, double>    RNGType;
   double volume;
   
   if (strcmp(vol_method,"sequence_of_balls")==0){ 
      if (strcmp(walk_method,"uniform_ball")==0){
         volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"CDHR")==0){
         volume = volume_sequence_of_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"RDHR")==0){
         volume = volume_sequence_of_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   }
   else if (strcmp(vol_method,"cooling_gaussian")==0){
      if (strcmp(walk_method,"gaussian_ball")==0){
         volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_CDHR")==0){
         volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_RDHR")==0){
         volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   } else if (strcmp(vol_method,"cooling_balls")==0){
       if (strcmp(walk_method,"uniform_ball")==0){
         volume = volume_cooling_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"CDHR")==0){
         volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"RDHR")==0){
         volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"billiard")==0){
         volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, epsilon, walk_len);
       }
   }
   return volume;
}
//////////          This is where the "COMPUTE_VOLUME" class ends  

//////////          This is where the "GENERATE_SAMPLES" class starts

//double HPolytopeCPP::generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, bool cdhr, bool rdhr, bool gaussian, bool set_L, bool billiard, bool ball_walk, double a, double L){
//double HPolytopeCPP::generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, char* boundary, char* cdhr, char* rdhr, char* gaussian, char* set_L, char* billiard, char* ball_walk, double a, double L){   
double HPolytopeCPP::generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, int boundary, int cdhr, int rdhr, int gaussian, int set_L, int billiard, int ball_walk, double a, double L){
   
   
   cout<<"Hello friend\n";
   
   if (boundary == 1) {
      cout<<"Friend your are doing well \n";
   } else {
      "What is wrong with you friend? \n";
   }
   
   double* samples;
   RNGType rng(HP.dimension());
   std::list<Point> rand_points;
   
   //Point default_starting_point = HP.ComputeInnerBall().first;
   Point starting_point = HP.ComputeInnerBall().first;
   
      std::cout<<"starting point with getCoefficients = "<<starting_point.getCoefficients()<<std::endl;
      std::cout<<"A = "<<HP.get_mat()<<std::endl;
      std::cout<<"b = "<<HP.get_vec()<<std::endl;
      std::cout<<"dimension = "<<HP.dimension()<<std::endl;
     
      std::cout<<"walk_len = "<<walk_len<<std::endl;
      std::cout<<"number_of_points = "<<number_of_points<<std::endl;
      //std::cout<<"starting_point = "<<starting_point<<std::endl;
      std::cout<<"number_of_points_to_burn = "<<number_of_points_to_burn<<std::endl;
      std::cout<<"walk_len = "<<walk_len<<std::endl; 
      
      std::cout<<"a = "<<a<<std::endl;
      std::cout<<"L = "<<L<<std::endl;   
   
   
   if (boundary == 1) {      
      if (cdhr == 1) {    
         uniform_sampling_boundary<BCDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
         } else {
            uniform_sampling_boundary<BRDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
         }
   } else if (cdhr == 1) {
      if (gaussian == 1) {
         gaussian_sampling<GaussianCDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<CDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      }
   } else if (rdhr == 1){
      if (gaussian == 1) {
         gaussian_sampling<GaussianRDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<RDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      }
   } else if (billiard == 1) {
      if (set_L == 1) {
         BilliardWalk WalkType(L);
         uniform_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<BilliardWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      }
   } else {
      if (set_L == 1) {
         if (gaussian == 1) {
            GaussianBallWalk WalkType(L);
            gaussian_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
            } else {
               BallWalk WalkType(L);
               uniform_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
            }
        } else {
            if (gaussian == 1) {
               gaussian_sampling<GaussianBallWalk>(rand_points, HP, rng, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
            } else {
               uniform_sampling<BallWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
            }
        }
   }
   
   std::cout<<"sampling completed"<<std::endl;

// The following block of code should NOT be removed!
   auto n_si=0;
   for (auto it_s = rand_points.begin(); it_s != rand_points.end(); it_s++){
      for (auto i = 0; i != it_s->dimension(); i++){
         samples[n_si++] = (*it_s)[i];
      }
   }
   
}
