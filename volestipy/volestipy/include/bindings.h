#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H

////// I think we could remove those...
//#include "Eigen/Eigen"
//#include "cartesian_kernel.h"

// from SOB volume - exactly the same for CG and CB methods
#include <fstream>
#include <iostream>
#include "random_walks.hpp"
#include "random.hpp"
#include "misc.h"
#include "known_polytope_generators.h"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

//from generate_samples, some extra headers not already included
#include <chrono>
#include "sampling/sampling.hpp"

// for rounding
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"
#include "preprocess/get_full_dimensional_polytope.hpp"    // will be used after completing the rounding



class HPolytopeCPP{
   public:

      typedef double NT;
      typedef Cartesian<NT>    Kernel;
      typedef typename Kernel::Point    Point;
      typedef HPolytope<Point> Hpolytope;
      typedef typename Hpolytope::MT    MT;
      typedef typename Hpolytope::VT    VT;
      typedef BoostRandomNumberGenerator<boost::mt19937, double>    RNGType;
      
      // regarding the rounding step
      typedef std::tuple<MT, VT, NT>    round_result;
      
      
      

// The class and its main specs
      HPolytopeCPP();
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables);

      Hpolytope HP;
      std::pair<Point,NT> CheBall;
      ~HPolytopeCPP();


// The functions of the HPoltopeCPP class

// the compute_volume() function
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);

// the generate_samples() function
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, \
       bool cdhr, bool rdhr, bool gaussian, bool set_L, bool billiard, bool ball_walk, double a, double L,  double* samples);

// the rounding() function
      void rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_value);

};




#endif
