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

//////////          This is where the "COMPUTE_VOLUME" class starts  

double HPolytopeCPP::compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed){
  
   typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
   double volume;
   
// strcmp returns a lexical difference (short-circuit serial byte comparator) of the two strings you have given as parameters. 0 means that both strings are equal  
// we are about to have 3 methods for computing the volume () and a number of random walks for each of those
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




////////          This is where the "GENERATE_SAMPLES" class starts  


//void HPolytopeCPP::generate_samples(double const& starting_point, unsigned int const& walk_len, unsigned int const& number_of_points, unsigned int const& number_of_points_to_burn ,
//                            bool const& boundary, bool const& cdhr, bool const& rdhr, bool const& gaussian, bool const& set_L, bool const& billiard, bool const& ball_walk,
//                            double const& a, double const& L){
//   
//   typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
//   typedef PointList &randPoints
//      
//   boost::random::uniform_real_distribution<>(urdist);
//   boost::random::uniform_real_distribution<> urdist1(-1,1);
//   std::list <Point> randPoints;
//   
//   
//   if (boundary == true) {
//      if (cdhr == true) {
//         uniform_sampling_boundary <BCDHRWalk>(randPoints, HP, RNGType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//         } else {
//            uniform_sampling_boundary <BRDHRWalk>(randPoints, HP, RNGType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//         }
//   } else if (cdhr == true) {
//      if (gaussian == true) {
//         gaussian_sampling<GaussianCDHRWalk>(randPoints, HP, RNGType, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
//      } else {
//         uniform_sampling<CDHRWalk>(randPoints, HP, RNGType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//      }
//   } else if (rdhr == true){
//      if (gaussian == true) {
//         gaussian_sampling<GaussianRDHRWalk>(randPoints, HP, RNGType, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
//      } else {
//         uniform_sampling<RDHRWalk>(randPoints, HP, RNGType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//      }
//   } else if (billiard == true) {
//      if (set_L == true) {
//         BilliardWalk WalkType(L);
//         uniform_sampling(randPoints, HP, RNGType, WalkType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//      } else {
//         uniform_sampling<BilliardWalk>(randPoints, HP, RNGType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//      }
//   } else {
//      if (set_L == true) {
//         if (gaussian == true) {
//            
//            GaussianBallWalk WalkType(L);
//            gaussian_sampling(randPoints, HP, RNGType, WalkType, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
//            
//            } else {
//               BallWalk WalkType(L);
//               uniform_sampling(randPoints, HP, RNGType, WalkType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//            }
//            
//        } else {
//            if (gaussian == true) {
//               gaussian_sampling<GaussianBallWalk>(randPoints, HP, RNGType, walk_len, number_of_points, a, starting_point, number_of_points_to_burn);
//            } else {
//               uniform_sampling<BallWalk>(randPoints, HP, RNGType, walk_len, number_of_points, starting_point, number_of_points_to_burn);
//            }
//        }
//   }
//
//// The following block of code should NOT be removed!
//   auto n_si=0;
//   for (auto it_s = randPoints.begin(); it_s != randPoints.end(); it_s++){
//      for (auto i = 0; i != it_s->dimension(); i++){
//         samples[n_si++] = (*it_s)[i];
//      }
//   }
//}
