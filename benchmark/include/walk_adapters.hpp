#pragma once

#include "core_types.hpp"
#include "benchmark_utils.hpp"
#include "walk_parameters.hpp"
#include "random_walks/random_walks.hpp" 

#include "sampling/sampling.hpp"
/*
 * -----------------
 * This file provides a uniform interface for all random walk implementations
 * used.
 *
 * Different walk algorithms expose different constructor signatures
 * and apply() methods. For example, most uniform samplers can be constructed
 * as Walk(P, p, rng), while Gaussian samplers require additional 
 * parameters (e.g. parameter a_i) and different apply() signatures.
 *
 * To avoid having special cases throughout the benchmark code, this file
 * introduces the WalkAdapter abstraction. Each adapter exposes a common
 * apply_batch() function that:
 *
 *   1. Constructs the corresponding walk object.
 *   2. Executes a specified number of walk steps.
 *   3. Collects generated sample points into a batch.
 *   4. Enforces benchmark time limits.
 * 
 * The overall philosophy is that every random walk method can be executed
 * through the same adapter interface, allowing the rest of the framework
 * (benchmark runners, statistics collection, batching logic, etc.) to remain
 * completely independent of the underlying sampling algorithm.
 * Check walk_run.hpp for the universal run function.
 * MACROS are used in an attempt to shorten the size of this file since most of the code is the same anyway
*/

// Uniform Walk Type definitions
typedef BallWalk::template Walk<HPOLYTOPE, RNGType> BallWalkType;
typedef BilliardWalk::template Walk<HPOLYTOPE, RNGType> BilliardWalkType;
typedef AcceleratedBilliardWalk::template Walk<HPOLYTOPE, RNGType> AcceleratedBilliardWalkType;
typedef SparseBilliardWalk::template Walk<HPOLYTOPE, RNGType> SparseBilliardWalkType;
typedef CDHRWalk::template Walk<HPOLYTOPE, RNGType> CDHRWalkType;
typedef RDHRWalk::template Walk<HPOLYTOPE, RNGType> RDHRWalkType;
typedef DikinWalk::template Walk<HPOLYTOPE, RNGType> DikinWalkType;
typedef JohnWalk::template Walk<HPOLYTOPE, RNGType> JohnWalkType;
typedef VaidyaWalk::template Walk<HPOLYTOPE, RNGType> VaidyaWalkType;

typedef GaussianBallWalk::template Walk<HPOLYTOPE, RNGType> GaussianBallWalkType;
typedef GaussianCDHRWalk::template Walk<HPOLYTOPE, RNGType> GaussianCDHRWalkType;

typedef ShakeAndBakeWalk::template Walk<HPOLYTOPE, RNGType> ShakeAndBakeWalkType;
typedef BilliardShakeAndBakeWalk::template Walk<HPOLYTOPE, RNGType> BilliardSBWalkType;

typedef BCDHRWalk::template Walk<HPOLYTOPE, RNGType> BCDHRWalkType;
typedef BRDHRWalk::template Walk<HPOLYTOPE, RNGType> BRDHRWalkType;


// For Uniform Walks
template <typename WalkType>
struct WalkAdapter {
    
    static constexpr bool supports_chunking = true;
    
    // Initialization 
    static WalkType init(HPOLYTOPE& P, Point& p, const BenchmarkConfig& config, RNGType& rng) {
        return WalkType(P, p, rng);
    }
    
    // apply batch
    static void apply_batch(WalkType& walk, HPOLYTOPE& P, Point& p, unsigned int batch_size, 
                            unsigned int walk_len, std::vector<Point>& batchPoints, 
                            const BenchmarkConfig& config, RNGType& rng, Timer& walk_timer) 
    {
        for (unsigned int i = 0; i < batch_size; ++i) {
            walk.apply(P, p, walk_len, rng); 
            batchPoints.push_back(p);

            if (i % 50 == 0 && walk_timer.get_total_time() > config.time_limit_sec) {
                break; 
            }
        }
    }
};

// For GAUSSIAN WALKs
#define REGISTER_GAUSSIAN_ADAPTER(WALK_TYPE, JSON_NAME) \
template <> \
struct WalkAdapter<WALK_TYPE> { \
    \
    static constexpr bool supports_chunking = true; \
    \
    /* Initialization */ \
    static WALK_TYPE init(HPOLYTOPE& P, Point& p, const BenchmarkConfig& config, RNGType& rng) { \
        double a_i = 1.0; \
        auto walk_iter = config.walk_settings.find(JSON_NAME); \
        if (walk_iter != config.walk_settings.end()) { \
            a_i = walk_iter->second.a_i_param; \
        } \
        return WALK_TYPE(P, p, a_i, rng); \
    } \
    \
    /* aplly batch */ \
    static void apply_batch(WALK_TYPE& walk, HPOLYTOPE& P, Point& p, unsigned int batch_size, \
                            unsigned int walk_len, std::vector<Point>& batchPoints, \
                            const BenchmarkConfig& config, RNGType& rng, \
                            Timer& walk_timer) \
    { \
        /* We fetch a_i to pass into the apply method */ \
        double a_i = 1.0; \
        auto walk_iter = config.walk_settings.find(JSON_NAME); \
        if (walk_iter != config.walk_settings.end()) { \
            a_i = walk_iter->second.a_i_param; \
        } \
        \
        for (unsigned int i = 0; i < batch_size; ++i) { \
            walk.apply(P, p, a_i, walk_len, rng); \
            batchPoints.push_back(p); \
            \
            /* The Timeout Check */ \
            if (i % 50 == 0 && walk_timer.get_total_time() > config.time_limit_sec) { \
                break; \
            } \
        } \
    } \
};

REGISTER_GAUSSIAN_ADAPTER(GaussianBallWalkType, "GaussianBallWalk")
REGISTER_GAUSSIAN_ADAPTER(GaussianCDHRWalkType, "GaussianCDHRWalk")

// FOR Billiard Shake-and-Bake walks
// For these walks we need a point on a facet to start with.
// To that end, we use the billiard logic to shoot a ray from our interior initial point and check where we hit a facet
template <>
struct WalkAdapter<BilliardSBWalkType> {

    static constexpr bool supports_chunking = true;

    // Initialization
    static BilliardSBWalkType init(HPOLYTOPE& P, Point& p, const BenchmarkConfig& config, RNGType& rng) {
        int nr = 10; // Hardcoded number of reflections 

        unsigned int n = P.dimension();
        Point v = GetDirection<Point>::apply(n, rng);

        // Temporary structures required by line_positive_intersect
        typename Point::Coeff lambdas(P.num_of_hyperplanes());
        typename Point::Coeff Av(P.num_of_hyperplanes());
        lambdas.setZero();
        Av.setZero();

        // Find intersection 
        std::pair<typename Point::FT, int> pbpair = P.line_positive_intersect(p, v, lambdas, Av);

        // Permanently move the starting point to the boundary
        p += (pbpair.first * v);
        int initial_facet = pbpair.second;

        // Construct and return the walk
        return BilliardSBWalkType(P, p, rng, initial_facet, nr);
    }

    // apply batch
    static void apply_batch(
        BilliardSBWalkType& walk,
        HPOLYTOPE& P,
        Point& p,
        unsigned int batch_size,
        unsigned int walk_len,
        std::vector<Point>& batchPoints,
        const BenchmarkConfig& config,
        RNGType& rng,
        Timer& walk_timer)
    {
        for (unsigned int i = 0; i < batch_size; ++i) {
            
            // Advance the already-initialized walk
            walk.apply(P, walk_len, rng);
            
            p = walk.getCurrentPoint();
            batchPoints.push_back(p);

            // Timeout Check
            if (i % 50 == 0 && walk_timer.get_total_time() > config.time_limit_sec) {
                break;
            }
        }
    }
};

// FOR Shake-and-Bake walks
template <>
struct WalkAdapter<ShakeAndBakeWalkType> {

    static constexpr bool supports_chunking = true;

    // Initialization. Find the boundary point
    static ShakeAndBakeWalkType init(HPOLYTOPE& P, Point& p, const BenchmarkConfig& config, RNGType& rng) {
        int initial_facet = 0;
        auto b = P.get_vec();
        bool on_boundary = false;

        // Check if the current point 'p' is already on a facet 
        for (int i = 0; i < P.num_of_hyperplanes(); ++i) {
            if (std::abs(P.get_row(i).dot(p.getCoefficients()) - b(i)) < 1e-7) {
                initial_facet = i;
                on_boundary = true;
                break;
            }
        }

        // If it's an interior point, project it to the boundary 
        if (!on_boundary) {
            typename Point::Coeff v_vec = Point::Coeff::Zero(P.dimension());
            v_vec(0) = 1.0; 

            double min_lambda = std::numeric_limits<double>::max();

            for (int i = 0; i < P.num_of_hyperplanes(); ++i) {
                double v_dot_a = P.get_row(i).dot(v_vec);
                if (v_dot_a > 1e-10) {
                    double dist = b(i) - P.get_row(i).dot(p.getCoefficients());
                    double lam = dist / v_dot_a;

                    if (lam > 0 && lam < min_lambda) {
                        min_lambda = lam;
                        initial_facet = i;
                    }
                }
            }

            // Move the initial point permanently to the collision spot
            typename Point::Coeff new_coords = p.getCoefficients() + (min_lambda * v_vec);
            p = Point(new_coords); 
        }

        // Construct and return the walk
        return ShakeAndBakeWalkType(P, p, initial_facet, rng);
    }

    // Update the walk
    static void apply_batch(ShakeAndBakeWalkType& walk, HPOLYTOPE& P, Point& p, 
                            unsigned int batch_size, unsigned int walk_len, 
                            std::vector<Point>& batchPoints, const BenchmarkConfig& config, 
                            RNGType& rng, Timer& walk_timer)
    {
        for (unsigned int i = 0; i < batch_size; ++i) {
           
            walk.apply(P, walk_len, rng);

            p = walk.getCurrentPoint();
            batchPoints.push_back(p);

            if (i % 50 == 0 && walk_timer.get_total_time() > config.time_limit_sec) {
                break;
            }
        }
    }
};


// FOR BOUNDARY HIT-AND-RUN WALKS (BCDHR, BRDHR)
#define REGISTER_BOUNDARY_HR_ADAPTER(WALK_TYPE, JSON_NAME) \
template <> \
struct WalkAdapter<WALK_TYPE> { \
    \
    static constexpr bool supports_chunking = true;  \
    \
    /* Initialization */ \
    static WALK_TYPE init(HPOLYTOPE& P, Point& p, const BenchmarkConfig& config, RNGType& rng) { \
        /* Constructor with 3 arguments */ \
        return WALK_TYPE(P, p, rng); \
    } \
    \
    /* Apply batch */ \
    static void apply_batch(WALK_TYPE& walk, HPOLYTOPE& P, Point& p, unsigned int batch_size, \
                            unsigned int walk_len, std::vector<Point>& batchPoints, \
                            const BenchmarkConfig& config, RNGType& rng, \
                            Timer& walk_timer) \
    { \
        /* Dummy points to catch the boundary chord endpoints for this batch */ \
        Point chord_p1 = p; \
        Point chord_p2 = p; \
        \
        for (unsigned int i = 0; i < batch_size; ++i) { \
            /* Advance the preserved walk object */ \
            walk.apply(P, chord_p1, chord_p2, walk_len, rng); \
            \
            /* Extract the updated internal point using getter */ \
            p = walk.getCurrentPoint(); \
            batchPoints.push_back(p); \
            \
            /* Timeout Check */ \
            if (i % 50 == 0 && walk_timer.get_total_time() > config.time_limit_sec) { \
                break; \
            } \
        } \
    } \
};

REGISTER_BOUNDARY_HR_ADAPTER(BCDHRWalkType, "BCDHRWalk")
REGISTER_BOUNDARY_HR_ADAPTER(BRDHRWalkType, "BRDHRWalk")

// Riemannian Hamiltonian
// This method works a bit differently so we run a big chunk of points at once.
template <>
struct WalkAdapter<CRHMCWalk> {
    
    static constexpr bool supports_chunking = false;

    // Dummy initialization so that we satisfy the generic template caller
    static int init(HPOLYTOPE& P, Point& p, const BenchmarkConfig& config, RNGType& rng) {
        return 0; 
    }

    // apply_batch takes the dummy integer, which we safely ignore
    static void apply_batch(int& dummy_state, HPOLYTOPE& P, Point& p, unsigned int batch_size, 
                            unsigned int walk_len, std::vector<Point>& batchPoints, 
                            const BenchmarkConfig& config, RNGType& rng, Timer& walk_timer) 
    {
        using Func = ZeroScalarFunctor<Point>;
        using Grad = ZeroFunctor<Point>;
        using Hess = ZeroFunctor<Point>;

        Func f;
        Grad g;
        Hess h;

        std::list<Point> temp_list;
        int n_burns = 1000; 

        // execute_crhmc handles the massive chunk and burn-in all at once
        execute_crhmc<HPOLYTOPE, RNGType, std::list<Point>, Grad, Func, Hess, CRHMCWalk, 1>(
            P, rng, temp_list, 1, batch_size, n_burns, &g, &f, &h
        );

        batchPoints.insert(batchPoints.end(), temp_list.begin(), temp_list.end());

        if (!temp_list.empty()) {
            p = temp_list.back();
        }
    }
};

