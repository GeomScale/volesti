#ifndef VOLESTIBINDINGS_H
#define VOLESTIBINDINGS_H


#include "Eigen/Eigen"


// altered by Haris at 2020_06_08
//#include "volume.h"
//#include "sampling.hpp"
#include "cartesian_kernel.h"



class HPolytopeCPP{
   public:

// ATTENTION! Only one of the following should remain. I need to talk to Tolis about this.
      //template <typename NT>;
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
      //double compute_volume(char* method, int walk_len, double epsilon, uint seed);
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, uint seed);
      void generate_samples(int walk_len, int n_samples, double* samples, uint seed);

};


#endif
