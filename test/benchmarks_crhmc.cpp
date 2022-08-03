// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// VolEsti example

#include "Eigen/Eigen"
#include "generators/known_polytope_generators.h"
#include <fstream>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "convex_bodies/hpolytope.h"
#include "preprocess/crhmc/crhmcProblem.h"
#include "preprocess/crhmc/crhmc_input.h"
#include "cartesian_geom/cartesian_kernel.h"
void benchmark(std::string fileName){
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  read_pointset(inp, Pin);
  inp.close();
  PolytopeType HP(Pin);
  int d = HP.dimension();
  Input input = Input(d);
  input.Aineq = HP.get_mat();
  input.bineq = HP.get_vec();
  double tstart;
  std::cout << "CRHMC polytope preparation for " <<fileName  << std::endl;
  tstart = (double)clock()/(double)CLOCKS_PER_SEC;
  CrhmcProblem P = CrhmcProblem(input);
  std::cout << "Preparation completed in time, ";
  std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart;
  std::cout<<" the resulting matrix has "<<P.Asp.nonZeros()<<" nonZeros"<< std::endl;
}
int main()
{   typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    using CrhmcProblem=crhmcProblem<Point>;
    using Input=crhmc_input<MT, NT>;

/*
    std::cout << "CRHMC polytope preparation" << std::endl << std::endl;

    Hpolytope HP = generate_cube<Hpolytope>(1000, false);
    int d = HP.dimension();
    Input input = Input(d);
    input.Aineq = HP.get_mat();
    input.bineq = HP.get_vec();
    double tstart;

    std::cout << "CRHMC polytope preparation (cube-1000)" << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    CrhmcProblem P = CrhmcProblem(input);

    std::cout << "Preparation completed in time, ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
*/
benchmark("../../test/metabolic_full_dim/polytope_e_coli.ine");
benchmark("../../test/netlib/afiro.ine");

return 0;
}
