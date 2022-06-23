// Using LDLT decomposition: more numerically stable for singular matrices
bool isPosSemidefinite(MT A){
    Eigen::LDLT<MT> A_ldlt(A);
    if (A_ldlt.info() != Eigen::NumericalIssue && A_ldlt.isPositive())
        return true;
    return false;
}

#include "matrix_operations/EigenvaluesProblems.h"
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/spectrahedra/spectrahedron.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
// typedef std::vector<double> VT;

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> spectrahedron;
typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
// typedef BilliardWalk::template Walk<spectrahedron, RNGType> BilliardWalkType;
// typedef AcceleratedBilliardWalk::template Walk<spectrahedron, RNGType> AcceleratedBilliardWalkType;

// Auxiliary geometric functions
Point getDirection(unsigned int const& dim, RNGType &rng, bool normalize=true){
    double normal = 0.;
    Point p(dim);
    double* data = p.pointerToData();

    for (unsigned int i=0; i<dim; ++i){
        *data = rng.sample_ndist();
        normal += *data * *data;
        data++;
    }

    normal = 1./std::sqrt(normal);
    if (normalize) p *= normal;
    return p;
}

// Sample uniformly distributed correlation matrices with Billiard Walk. 
// Report the n that the sampler takes more than 3 minutes to sample 1000 matrices

// This function is taken and simplified from uniform_billiard_walk.hpp

Point BilliardWalkSpectra(spectrahedron &P, Point& q, unsigned int const& walk_length, unsigned int nreflex, RNGType &rng, double const _Len){
    unsigned int k = P.dimension();
    double L, tau;
    Point p = q, v;
    std::pair<double, int> pbpair;
    for (unsigned int j=0; j<walk_length; ++j){
        L = rng.sample_urdist() * _Len;
        v = getDirection(k, rng);
        Point p0 = p;
        int it = 0;
        while (it < nreflex)
        {
            pbpair = intersection(P, p, v, k);
            tau = pbpair.first;
            if (L <= tau) {
                p += (L * v);
                break;
            }
            tau = 0.995 * tau; // 0.995: to approximate boundary points?
            p += tau * v; // A point (almost) on the boundary
            L -= tau;
            reflection(P, p, v, pbpair.second);
            it++;
        }
        if (it == nreflex){
            p = p0;
        }
    }
    return p;
}

std::vector<Point> uniform_correl_matrix_sampling(unsigned int n, unsigned int num_points, unsigned int walk_len=10, unsigned int nreflex = 10){
    int k = n*(n-1)/2; // Dimension: k = n(n-1)/2
    RNGType rng(k);
    Point p(k); // Initial interior point (origin: identity matrix)
    std::vector<Point> randPoints;

    // Create the spectrahedron
    std::vector<MT> lmi_mat = myLMIGenerator(n);
    LMI<double, MT, VT> lmi(lmi_mat);
    spectrahedron spectra(lmi);
    spectra.set_interior_point(p);
    std::pair<Point, double> inner_ball = spectra.ComputeInnerBall();
    double diameter = 6 * k * inner_ball.second;
    
    for (unsigned int i = 0; i < num_points; ++i){
        p = BilliardWalkSpectra(spectra, p, walk_len, nreflex, rng, diameter);
        randPoints.push_back(p);
    }
    return randPoints;
}

// Test 3:

int Test3(unsigned int num_points=1000, unsigned int walk_len=10, unsigned int nreflex = 10){
    int n = 3;
    double time;
    while(true){

        auto start = std::chrono::steady_clock::now();

        uniform_correl_matrix_sampling(n, num_points, walk_len, nreflex);

        auto end = std::chrono::steady_clock::now();
        time = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
        std::cout << time << std::endl;
        if(time > 3) return n;
        ++n;
    }
}

template <typename NT, typename WalkType>
MT uniform_sampling_correlation_matrices(unsigned int n){
    CorreSpectra<Point> P(n);

    MT samples
    switch(WalkType){
        case BilliardWalk:

        case UniformHamiltonianMonteCarloWalk:

    }
    return samples;
}

template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename Point
>
void uniform_sampling(PointList &randPoints,
                   Polytope &P,
                   RandomNumberGenerator &rng,
                   WalkTypePolicy &WalkType,
                   const unsigned int &walk_len,
                   const unsigned int &rnum,
                   const Point &starting_point,
                   unsigned int const& nburns)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;
    typedef RandomPointGenerator<walk> RandomPointGenerator;

    Point p = starting_point;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, rnum, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}

template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename NT,
        typename Point
        >
void gaussian_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       WalkTypePolicy &WalkType,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef GaussianRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, a, rnum, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}

template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename NT,
        typename Point,
        typename NegativeGradientFunctor,
        typename NegativeLogprobFunctor,
        typename Solver
        >
void logconcave_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const Point &starting_point,
                       unsigned int const& nburns,
                       NegativeGradientFunctor &F,
                       NegativeLogprobFunctor &f)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Point,
                    Polytope,
                    RandomNumberGenerator,
                    NegativeGradientFunctor,
                    NegativeLogprobFunctor,
                    Solver
            > walk;

    typedef typename WalkTypePolicy::template parameters
            <
                    NT,
                    NegativeGradientFunctor
            > walk_params;

    // Initialize random walk parameters
    unsigned int dim = starting_point.dimension();
    walk_params params(F, dim);

    if (F.params.eta > 0) {
        params.eta = F.params.eta;
    }

    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    walk logconcave_walk(&P, p, F, f, params);

    typedef LogconcaveRandomPointGenerator<walk> RandomPointGenerator;

    RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                push_back_policy, rng, F, f, params, logconcave_walk);

    logconcave_walk.disable_adaptive();
    randPoints.clear();

    RandomPointGenerator::apply(P, p, rnum, walk_len, randPoints,
                                push_back_policy, rng, F, f, params, logconcave_walk);

}

template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename NT,
        typename Point
        >
void exponential_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       WalkTypePolicy &WalkType,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const Point &c,
                       const NT &a,
                       const NT &eta,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef ExponentialRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, c, a, eta, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, c, a, eta, rnum, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}