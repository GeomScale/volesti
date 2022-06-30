#include <fstream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "known_polytope_generators.h"
#include "sampling/sampling.hpp"

#include "diagnostics/univariate_psrf.hpp"

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

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> spectrahedron;
typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
// typedef BilliardWalk::template Walk<spectrahedron, RNGType> BilliardWalkType;
// typedef AcceleratedBilliardWalk::template Walk<spectrahedron, RNGType> AcceleratedBilliardWalkType;

// Test 2: LMI Generator for Test 3

std::vector<MT> myLMIGenerator(int n){
    int i, j, l, k = n*(n-1)/2+1;
    std::vector<MT> list_Mat;
    MT A;
    list_Mat.push_back(MT::Identity(n, n));
    for(i = 0; i < n; i++){
        for(j = i+1; j < n; j++){
            A = MT::Zero(n, n);
            A(i,j) = 1;
            A(j,i) = 1;
            list_Mat.push_back(A);
        }
    }
    return list_Mat;
}

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

bool membership(spectrahedron &spectra, const VT &xvector, const unsigned int n, const unsigned int k){
    // std::vector<double>::iterator it = xvector.begin();
    for(int i = 0; i < k; ++i)
        if((xvector(i) > 1) || (xvector(i) < -1)) return false;
    MT A = rebuildMatrix(xvector, n);
    if(isPosSemidefinite(A)) return true;
    return false;
}

std::pair<double, int> intersection(spectrahedron &P, const Point &x, const Point &v, const unsigned int k){
    double tau, tmp;
    int j = 0;
    if(v[0] > 0){
        tau = (1-x[0])/v[0];   
    }else{
        tau = -(1 + x[0])/v[0];
    }
    for(int i = 1; i < k; ++i){
        if(v[i] > 0){
            tmp = (1 - x[i])/v[i];
        }else{
            tmp = -(1 + x[i])/v[i];
        }
        if(tau > tmp){
            tau = tmp;
            j = i;
        }
    }
    tmp = P.positiveLinearIntersection(x.getCoefficients(), v.getCoefficients());
    if(tau > tmp){
        tau = tmp;
        j = -1;
    }
    std::pair<double, int> res(tau,j);
    return res;
}

void reflection(spectrahedron P, Point &p, Point &v, const int flag){
    if(flag != -1){
        v.set_coord(flag, - v.getCoefficients()(flag));
        return;
    }
    P.compute_reflection(v, p, flag);
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

template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT sample_correlation_matrices(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;
    NT a = 1;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();

    int d = n*(n-1)/2;
    
    RNGType rng(d);
    Point StartingPoint(d);

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

    std::list<Point> randPoints;

    gaussian_sampling<WalkType>(randPoints, P, rng, walkL, numpoints, a,
                               StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

std::vector<Point> uniform_correl_matrix_sampling(unsigned int n, unsigned int num_points, unsigned int walk_len=10, unsigned int nreflex = 10){
    int k = n*(n-1)/2; // Dimension: k = n(n-1)/2
    RNGType rng(k);
    Point p(k); // Initial interior point (origin: identity matrix)
    std::vector<Point> randPoints;

    correlation_matrix matrix(lmi);
    spectra.set_interior_point(p);
    std::pair<Point, double> inner_ball = spectra.ComputeInnerBall();
    double diameter = 6 * k * inner_ball.second;
    
    for (unsigned int i = 0; i < num_points; ++i){
        p = BilliardWalkSpectra(spectra, p, walk_len, nreflex, rng, diameter);
        randPoints.push_back(p);
    }
    return randPoints;
}

NT radius = maxDouble;

        for (unsigned int i = 0; i < dimension(); ++i) {

            std::pair<NT, NT> min_max = coordinateIntersection(_inner_ball.first.getCoefficients(), i+1);

            if (min_max.first < radius) radius = min_max.first;
            if (-min_max.second < radius) radius = -min_max.second;
        }

        radius = radius / std::sqrt(NT(dimension()));
        _inner_ball.second = radius;

        return std::pair<PointType, NT>(_inner_ball.first, radius);

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
    typedef RandomPointGenerator<walk> RandomPointGenerator;

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

void exponential_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       WalkTypePolicy &WalkType,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const Point &c,
                       const NT &a,
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
        RandomPointGenerator::apply(P, p, c, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
                                    RandomPointGenerator::apply(P, p, c, a, eta, nburns, walk_len, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, c, a, rnum, walk_len, randPoints,
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