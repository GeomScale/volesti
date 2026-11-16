#pragma once

#include <list>
#include <stdexcept>

#include "core_types.hpp"
#include "benchmark_utils.hpp"
#include "walk_parameters.hpp"
#include "random_walks/random_walks.hpp"
#include "random_walks/shake_and_bake_walk.hpp"
#include "random_walks/billiard_shake_and_bake_walk.hpp" 

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
        for (unsigned int i = 0; i < batch_size / 2; ++i) { \
            /* Advance the preserved walk object */ \
            walk.apply(P, chord_p1, chord_p2, walk_len, rng); \
            \
            /* Extract the boundary points*/ \
            batchPoints.push_back(chord_p1); \
            batchPoints.push_back(chord_p2); \
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

// Riemannian Hamiltonian (CRHMC)
// Preprocessing + burn-in happen in init() while apply_batch() only advances the chain.
template <>
struct WalkAdapter<CRHMCWalk> {

    using Func   = ZeroScalarFunctor<Point>;
    using Grad   = ZeroFunctor<Point>;
    using Hess   = ZeroFunctor<Point>;
    using Input  = crhmc_input<typename HPOLYTOPE::MT, Point, Func, Grad, Hess>;
    using Problem = crhmc_problem<Point, Input>;
    using Solver = ImplicitMidpointODESolver<Point, NT, Problem, Grad, 1>;
    using Walk   = CRHMCWalk::template Walk<Point, Problem, RNGType, Grad, Func, Solver>;
    using Params = CRHMCWalk::template parameters<NT, Grad>;
    using Generator = CrhmcRandomPointGenerator<Walk>;

    struct State {
        Func f; Grad g; Hess h;
        std::unique_ptr<Input>   input;
        std::unique_ptr<Problem> problem;
        std::unique_ptr<Params>  params;
        std::unique_ptr<Point>   p;
        std::unique_ptr<Walk>    walk;
    };

    static std::unique_ptr<State> init(HPOLYTOPE& P, Point& /*p*/, const BenchmarkConfig& config, RNGType& rng) {
        auto s = std::make_unique<State>();

        // Same argument order as in crhmc_sampling: (P, f, F, h)
        s->input = std::make_unique<Input>(
            convert2crhmc_input<Input, HPOLYTOPE, Func, Grad, Hess>(P, s->f, s->g, s->h));

        s->problem = std::make_unique<Problem>(*s->input);
        if (s->problem->terminate) {
            // Failing loudly avoids an endless chunk loop with zero points generated
            throw std::runtime_error("CRHMC preprocessing terminated (problem.terminate == true)");
        }

        s->p = std::make_unique<Point>(s->problem->center);
        s->problem->options.simdLen = 1;

        s->params = std::make_unique<Params>(s->input->df, s->p->dimension(), s->problem->options);
        if (s->input->df.params.eta > 0) {
            s->params->eta = s->input->df.params.eta;
        }

        s->walk = std::make_unique<Walk>(*s->problem, *s->p, s->input->df, s->input->f, *s->params);

        int n_burns = 1000;
        std::list<Point> burn;
        PushBackWalkPolicy push_back_policy;
        Generator::apply(*s->problem, *s->p, n_burns, 1, burn, push_back_policy,
                         rng, s->g, s->f, *s->params, *s->walk);

        return s;
    }

    static void apply_batch(std::unique_ptr<State>& s, HPOLYTOPE& P, Point& p,
                            unsigned int batch_size, unsigned int walk_len,
                            std::vector<Point>& batchPoints, const BenchmarkConfig& config,
                            RNGType& rng, Timer& walk_timer)
    {
        // The generator can't be interrupted mid-call, so the time check happens
        // before each chunk
        if (walk_timer.get_total_time() > config.time_limit_sec) return;

        std::list<Point> out;
        PushBackWalkPolicy push_back_policy;

        // Continues from *s->p with the same walk object. Returned points are already
        // mapped back to the original space (raw_output = false).
        Generator::apply(*s->problem, *s->p, batch_size, walk_len, out, push_back_policy,
                         rng, s->g, s->f, *s->params, *s->walk, 1, false);

        batchPoints.insert(batchPoints.end(), out.begin(), out.end());
        if (!out.empty()) p = out.back();   // p is only for the caller, the walk uses *s->p
    }
};

