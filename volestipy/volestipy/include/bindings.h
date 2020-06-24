#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H

//// I think we could remove those...
#include "Eigen/Eigen"
#include "cartesian_kernel.h"

// from SOB volume - exactly the same for CG and CB methods 
#include <fstream>
#include <iostream>
#include "random_walks.hpp"
#include "random.hpp"
#include "misc.h"
#include "known_polytope_generators.h"
#include "random/uniform_int.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

//from generate_samples, some extra headers not already included
#include <chrono>
#include "sampling/sampling.hpp"



class HPolytopeCPP{
   public:

      typedef double NT;
      typedef Cartesian<NT>    Kernel;
      typedef typename Kernel::Point    Point;
      typedef HPolytope<Point> Hpolytope;
      typedef typename Hpolytope::MT    MT;
      typedef typename Hpolytope::VT    VT;      
      typedef BoostRandomNumberGenerator<boost::mt19937, double>    RNGType;
      
      
      HPolytopeCPP();
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables);

      Hpolytope HP;
      std::pair<Point,NT> CheBall;
      ~HPolytopeCPP();
      
      //Point default_starting_point = HP.ComputeInnerBall().first;

      
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);      
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, bool cdhr, bool rdhr, bool gaussian, bool set_L, bool billiard, bool ball_walk, double a, double L);
      //double* default_starting_point    --> include it on previous line
};


#endif
