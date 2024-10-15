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

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;

template <typename Polytope>
void mmcs_sampling(Polytope const& P)
{
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT S;
    int total_neff;

    mmcs(P, 100, S, total_neff);

    for (int i = 0; i < S.cols(); ++i) {
        Point p(S.rows());
        for (int j = 0; j < S.rows(); ++j) {
            p.set_coord(j, S(j,i));
        }
        if (P.is_in(p) == 0) {
            std::cout << "Sample point out of the polytope.";
        }
    }

    std::cerr << "sum of effective sample sizes: " << total_neff << std::endl;
    std::cerr << "multivariate PSRF: " <<  multivariate_psrf<NT, VT>(S) << std::endl;
    std::cerr << "maximum marginal PSRF: " <<  univariate_psrf<NT, VT>(S).maxCoeff() << std::endl;
}

int main()
{
    typedef boost::mt19937 PolyRNGType;

    int n = 50;
    Hpolytope P1 = random_hpoly<Hpolytope, PolyRNGType>(n, 4*n, 127); // we fix the example polytope, seed = 127
    mmcs_sampling(P1);

    Hpolytope P2 = generate_skinny_cube<Hpolytope>(10, false);
    mmcs_sampling(P2);

    return 0;
}

