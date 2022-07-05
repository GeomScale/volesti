#include "matrix_operations/EigenvaluesProblems.h"
#include <vector>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "misc.h"

#include "random_walks/random_walks.hpp"

#include "sampling/sampling.hpp"

#include "diagnostics/univariate_psrf.hpp"

#include "CorreSpectra.cpp"

template<typename PointType,
    typename PointList,
    typename RandomNumberGenerator, 
    typename WalkTypePolicy>
std::vector<Point> uniform_correlation_sampling(PointList &randPoints,
                   int n,
                   RandomNumberGenerator &rng,
                   WalkTypePolicy &WalkType,
                   const unsigned int &walk_len,
                   const unsigned int &num_points,
                   unsigned int const& nburns){
    typedef typename WalkTypePolicy::template Walk <CorreSpectra,
                                                    RandomNumberGenerator> walk;
    typedef RandomPointGenerator <walk> RandomPointGenerator;

    PushBackWalkPolicy push_back_policy;
    int d = n*(n-1)/2;
    RandomNumberGenerator rng(d);
    CorreSpectra<Point> P(n);
    Point p(d);
    P.set_interior_point(p);
    
    
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, num_points, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}