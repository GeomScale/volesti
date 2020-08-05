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

   int index = 0;
   for (int i = 0; i < n_hyperplanes; i++){
      b(i) = b_np[i];
      for (int j=0; j < n_variables; j++){
         A(i,j) = A_np[index];
         index++;
      }
   }

   HP.init(n_variables,A,b);
   CheBall = HP.ComputeInnerBall();
}
HPolytopeCPP::~HPolytopeCPP(){}


//////////          start of "compute_volume"          //////////
double HPolytopeCPP::compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed){
  
   double volume;

   if (strcmp(vol_method,"sequence_of_balls") == 0){
      if (strcmp(walk_method,"uniform_ball") == 0){
         volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"CDHR") == 0){
         volume = volume_sequence_of_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"RDHR") == 0){
         volume = volume_sequence_of_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   }
   else if (strcmp(vol_method,"cooling_gaussian") == 0){
      if (strcmp(walk_method,"gaussian_ball") == 0){
         volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_CDHR") == 0){
         volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, epsilon, walk_len);
      } else if (strcmp(walk_method,"gaussian_RDHR") == 0){
         volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(HP, epsilon, walk_len);
      }
   } else if (strcmp(vol_method,"cooling_balls") == 0){
       if (strcmp(walk_method,"uniform_ball") == 0){
         volume = volume_cooling_balls<BallWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"CDHR") == 0){
         volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"RDHR") == 0){
         volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, epsilon, walk_len);
       } else if (strcmp(walk_method,"billiard") == 0){
         volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, epsilon, walk_len);
       }
   }
   return volume;
}
//////////           end of "compute_volume()"            //////////


//////////         start of "generate_samples()"          //////////
double HPolytopeCPP::generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, bool cdhr, bool rdhr, bool gaussian, bool set_L, bool billiard, bool ball_walk, double a, double L, double* samples){

   RNGType rng(HP.dimension());
   std::list<Point> rand_points;

   //Point default_starting_point = HP.ComputeInnerBall().first;
   Point starting_point = HP.ComputeInnerBall().first;

   if (boundary == true) {
      if (cdhr == true) {
         uniform_sampling_boundary<BCDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
         } else {
            uniform_sampling_boundary<BRDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
         }
   } else if (cdhr == true) {
      if (gaussian == true) {
         gaussian_sampling<GaussianCDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<CDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      }
   } else if (rdhr == true){
      if (gaussian == true) {
         gaussian_sampling<GaussianRDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<RDHRWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      }
   } else if (billiard == true) {
      if (set_L == true) {
         BilliardWalk WalkType(L);
         uniform_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      } else {
         uniform_sampling<BilliardWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
      }
   } else {
      if (set_L == true) {
         if (gaussian == true) {
            GaussianBallWalk WalkType(L);
            gaussian_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
            } else {
               BallWalk WalkType(L);
               uniform_sampling(rand_points, HP, rng, WalkType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
            }
        } else {
            if (gaussian == true) {
               gaussian_sampling<GaussianBallWalk>(rand_points, HP, rng, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
            } else {
               uniform_sampling<BallWalk>(rand_points, HP, rng, walk_len, number_of_points, starting_point, number_of_points_to_burn);
            }
        }
   }
   std::cout<<"Sampling completed"<<std::endl;

// The following block of code allows us to parse the matrix with the points we are making
   auto n_si=0;
   for (auto it_s = rand_points.begin(); it_s != rand_points.end(); it_s++){
      for (auto i = 0; i != it_s->dimension(); i++){
         samples[n_si++] = (*it_s)[i];
      }
   }
}
//////////         end of "generate_samples()"          //////////


//////////         start of "rounding()"          //////////

double HPolytopeCPP::rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_val){
//  the initial function returns a tuple with the followings: (T, shift, std::abs(round_val))

   auto P(HP);
   RNGType rng(HP.dimension());
   CheBall = P.ComputeInnerBall();

   std::cout<<"CheBall Coefficient equals to = "<<P.ComputeInnerBall().first.getCoefficients()<<std::endl;
   std::cout<<"CheBall second equals to = "<<P.ComputeInnerBall().second<<std::endl;
   
   round_result round_res;
   //double round_val;
   
   int walk_len = 2;
   
   if (strcmp(rounding_method,"min_ellipsoid") == 0){
      round_res = min_sampling_covering_ellipsoid_rounding<AcceleratedBilliardWalk, MT, VT>(P, CheBall, walk_len, rng);      
   } else if (strcmp(rounding_method,"svd") == 0){      
      round_res = svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, CheBall, walk_len, rng);   
   } else if (strcmp(rounding_method, "max_ellipsoid") == 0){
       round_res = max_inscribed_ellipsoid_rounding<MT, VT>(P, CheBall); 
   }
      
   MT A_to_copy = P.get_mat();
   int new_n_hyperplanes = A_to_copy.rows();
   int new_n_variables = P.dimension();
   
   auto n_si = 0;
   for (int i = 0; i < new_n_hyperplanes; i++){
      for (int j = 0; j < new_n_variables; j++){
         new_A[n_si++] = A_to_copy(i,j);
      }
   }   
   
   VT b_to_copy = P.get_vec();
   for (int i=0; i < new_n_hyperplanes; i++){
      new_b[i] = b_to_copy[i];   
   }

   MT T_to_copy = get<0>(round_res);
   auto t_si = 0;
   for (int i = 0; i < new_n_hyperplanes; i++){
      for (int j = 0; j < new_n_variables; j++){
         T_matrix[t_si++] = T_to_copy(i,j);
      }
   }
   
   VT shift_to_copy = get<1>(round_res);
   for (int i = 0; i < new_n_hyperplanes; i++){
      shift[i] = shift_to_copy[i];
   }
   
   round_val = get<2>(round_res);

   std::cout<<"Till here my memory is good.\n"<<std::endl;
   std:cout<< round_val << endl;
   
   return 1.0;

}
