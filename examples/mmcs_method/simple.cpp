// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

#include "Eigen/Eigen"

#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "sampling/mmcs.hpp"
#include "generators/h_polytopes_generator.h"
#include "generators/known_polytope_generators.h"
#include "diagnostics/multivariate_psrf.hpp"
#include "diagnostics/univariate_psrf.hpp"
#include "diagnostics/ess_window_updater.hpp"


template <typename NT>
void run_main()
{
    typedef Cartesian<NT> Kernel;
    typedef boost::mt19937 PolyRNGType;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    {
        int n = 50;
        Hpolytope P = random_hpoly<Hpolytope, PolyRNGType>(n, 4*n, 127); // we fix the example polytope, seed = 127

        MT S;
        int total_neff;
        mmcs(P, 400, S, total_neff);

        std::cerr << "sum of effective sample sizes: " << total_neff << std::endl;
        std::cerr << "multivariate PSRF: " <<  multivariate_psrf<NT, VT>(S) << std::endl;
        std::cerr << "maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(S).maxCoeff() << std::endl;
    }
    {
        int n = 10;
        Hpolytope P = generate_skinny_cube<Hpolytope>(n, false);

        MT S;
        int total_neff;
        mmcs(P, 1000, S, total_neff);

        std::cerr << "sum of effective sample sizes: " << total_neff << std::endl;
        std::cerr << "multivariate PSRF: " <<  multivariate_psrf<NT, VT>(S) << std::endl;
        std::cerr << "maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(S).maxCoeff() << std::endl;
    }
}

int main() {
  run_main<double>();
  return 0;
}
