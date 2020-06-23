#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H

//// I think we could remove those...
#include "Eigen/Eigen"
#include "cartesian_kernel.h"


// from SOB volume - exactly the same for CG and CB methods 
#include <fstream>
#include <iostream>

#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include "known_polytope_generators.h"

//#include "doctest.h"
#include "misc.h"

//from generate_samples, some extra headers not already included
#include <chrono>
#include "sampling/sampling.hpp"



class HPolytopeCPP{
   public:

      typedef double NT;

      typedef Cartesian<NT>    Kernel;
      typedef typename Kernel::Point    Point;
      typedef boost::mt19937    RNGType;
      typedef HPolytope<Point> Hpolytope;

      typedef typename Hpolytope::MT    MT;
      typedef typename Hpolytope::VT    VT;

      HPolytopeCPP();
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables);

      Hpolytope HP;
      std::pair<Point,NT> CheBall;
      ~HPolytopeCPP();

      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);
      
      //void generate_samples(double const& starting_point, unsigned int const& walk_len, unsigned int const& number_of_points, unsigned int const& number_of_points_to_burn,
      //                      bool const& boundary, bool const& cdhr, bool const& rdhr, bool const& gaussian, bool const& set_L, bool const& billiard, bool const& ball_walk,
      //                      double const& a, double const& L)
      
};


#endif
