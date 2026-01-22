/*
 * sample_points.cpp
 * 
 * GNU Octave interface for Volesti library
 * Point sampling from H-polytopes and V-polytopes
 * 
 * This implements uniform sampling using CDHR, RDHR, and Ball Walk algorithms.
 * Architecture follows volume.cpp with zero-copy Eigen::Map pattern.
 */

#include <octave/oct.h>
#include <octave/parse.h>
#include <Eigen/Eigen>
#include <list>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "volume/sampling_policies.hpp"
#include "sampling/random_point_generators.hpp"
#include "sampling/sphere.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "generators/boost_random_number_generator.hpp"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;
typedef VPolytope<Point> Vpolytope;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

enum WalkType {
    CDHR = 0,
    RDHR = 1,
    BALL_WALK = 2
};

DEFUN_DLD(sample_points, args, nargout,
          "samples = sample_points(A, b, n [, walk_type, walk_length, nburns, verbose])\n\
          samples = sample_points(V, n [, walk_type, walk_length, nburns, verbose])\n\
\n\
Sample points uniformly from a polytope using MCMC.\n\
\n\
H-Polytope (Ax <= b):\n\
  A : m x n matrix of constraint coefficients\n\
  b : m x 1 vector of constraint bounds\n\
\n\
V-Polytope (convex hull of vertices):\n\
  V : m x n matrix where each row is a vertex\n\
\n\
Parameters (both types):\n\
  n : number of samples to generate\n\
  walk_type : (optional) 0=CDHR, 1=RDHR, 2=BallWalk (default: auto)\n\
  walk_length : (optional) steps per sample (default: 1)\n\
  nburns : (optional) burn-in iterations (default: 0)\n\
  verbose : (optional) show progress (default: true)\n\
\n\
Returns:\n\
  samples : d x n matrix, each column is a sample point\n\
\n\
Examples:\n\
  % H-polytope: 2D unit square\n\
  A = [1 0; -1 0; 0 1; 0 -1];\n\
  b = ones(4, 1);\n\
  samples = sample_points(A, b, 100);\n\
\n\
  % V-polytope: 3D simplex\n\
  V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];\n\
  samples = sample_points(V, 100);\n\
\n\
  % Custom walk parameters\n\
  samples = sample_points(A, b, 100, 0, 10, 50);  % CDHR, walk_length=10, nburns=50\n")
{
    if (args.length() < 2 || args.length() > 7)
    {
        error("sample_points: 2 to 7 arguments required");
        return octave_value_list();
    }

    // Determine polytope type based on arguments
    // H-polytope: (A, b, n, ...) where b is column vector
    // V-polytope: (V, n, ...) where n is scalar
    bool is_hpolytope = false;
    bool is_vpolytope = false;
    
    if (args.length() >= 2)
    {
        // Check if second argument is scalar (V-poly with n) or vector (H-poly with b)
        if (args(1).is_scalar_type())
        {
            is_vpolytope = true;
        }
        else
        {
            is_hpolytope = true;
        }
    }

    // Parse parameters (shifted index for V-polytope)
    int param_offset = is_hpolytope ? 3 : 2;  // After (A,b,n) or (V,n)
    
    unsigned int numpoints;
    int walk_type = -1;  // -1 = auto-select
    unsigned int walk_length = 1;
    unsigned int nburns = 0;
    bool verbose = true;
    
    // Get number of samples
    if (is_hpolytope)
    {
        if (args.length() < 3)
        {
            error("sample_points: H-polytope requires at least 3 arguments (A, b, n)");
            return octave_value_list();
        }
        numpoints = args(2).scalar_value();
    }
    else
    {
        numpoints = args(1).scalar_value();
    }
    
    if (numpoints <= 0)
    {
        error("sample_points: Number of samples must be positive");
        return octave_value_list();
    }
    
    // Parse optional parameters
    if (args.length() >= param_offset + 1)
    {
        walk_type = args(param_offset).scalar_value();
        if (walk_type < 0 || walk_type > 2)
        {
            error("sample_points: walk_type must be 0 (CDHR), 1 (RDHR), or 2 (BallWalk)");
            return octave_value_list();
        }
    }
    
    if (args.length() >= param_offset + 2)
    {
        double walk_length_dbl = args(param_offset + 1).scalar_value();
        if (walk_length_dbl < 1 || walk_length_dbl > 1e6)
        {
            error("sample_points: walk_length must be between 1 and 1e6");
            return octave_value_list();
        }
        walk_length = static_cast<unsigned int>(walk_length_dbl);
    }
    
    if (args.length() >= param_offset + 3)
    {
        double nburns_dbl = args(param_offset + 2).scalar_value();
        if (nburns_dbl < 0 || nburns_dbl > 1e6)
        {
            error("sample_points: nburns must be between 0 and 1e6");
            return octave_value_list();
        }
        nburns = static_cast<unsigned int>(nburns_dbl);
    }
    
    if (args.length() >= param_offset + 4)
    {
        verbose = args(param_offset + 3).bool_value();
    }

    // Main sampling logic
    Matrix result_matrix;
    
    if (is_hpolytope)
    {
        // H-Polytope sampling
        Matrix octave_A;
        ColumnVector octave_b;
        int m, n;
        
        try
        {
            octave_A = args(0).matrix_value();
            octave_b = args(1).column_vector_value();
            
            m = octave_A.rows();    // number of constraints
            n = octave_A.cols();    // dimension
        }
        catch (const std::exception& e)
        {
            error("sample_points: Invalid input types. Expected matrix A and column vector b. Error: %s", e.what());
            return octave_value_list();
        }
        
        if (octave_b.numel() != m)
        {
            error("sample_points: Dimensions of A and b do not match");
            return octave_value_list();
        }
        
        if (n < 1 || m < n + 1)
        {
            error("sample_points: Invalid polytope dimensions (need n >= 1 and m >= n+1, got n=%d, m=%d)", n, m);
            return octave_value_list();
        }
        
        // Auto-select walk type for H-polytope
        if (walk_type == -1)
        {
            walk_type = CDHR;  // Default for H-polytope
        }
        
        // Map Octave data to Eigen (zero-copy)
        Eigen::Map<MT> A_eigen(octave_A.fortran_vec(), m, n);
        Eigen::Map<VT> b_eigen(octave_b.fortran_vec(), m);
        
        Hpolytope P(n, A_eigen, b_eigen);
        
        if (verbose)
        {
            octave_stdout << "[Volesti C++] Sampling " << numpoints << " points from " 
                          << n << "D H-polytope (" << m << " constraints)..." << std::endl;
            const char* walk_names[] = {"CDHR", "RDHR", "BallWalk"};
            octave_stdout << "[Volesti C++] Walk: " << walk_names[walk_type]
                          << ", walk_length=" << walk_length 
                          << ", nburns=" << nburns << std::endl;
        }
        
        // Compute inner ball for starting point
        std::pair<Point, NT> inner_ball;
        try
        {
            inner_ball = P.ComputeInnerBall();
            if (inner_ball.second < 0.0)
            {
                error("sample_points: Unable to compute a feasible starting point");
                return octave_value_list();
            }
        }
        catch (const std::exception& e)
        {
            error("sample_points: Failed to compute inner ball: %s", e.what());
            return octave_value_list();
        }
        
        Point starting_point = inner_ball.first;
        
        // Create RNG
        RNGType rng(n);
        
        // Sample points using inline RandomPointGenerator pattern
        std::list<Point> rand_points;
        PushBackWalkPolicy push_back_policy;
        Point p = starting_point;
        
        try
        {
            switch (walk_type)
            {
            case CDHR:
                {
                    typedef typename CDHRWalk::template Walk<Hpolytope, RNGType> WalkType;
                    typedef RandomPointGenerator<WalkType> PointGenerator;
                    
                    if (nburns > 0) {
                        PointGenerator::apply(P, p, nburns, walk_length, rand_points,
                                            push_back_policy, rng);
                        rand_points.clear();
                    }
                    PointGenerator::apply(P, p, numpoints, walk_length, rand_points,
                                        push_back_policy, rng);
                }
                break;
            case RDHR:
                {
                    typedef typename RDHRWalk::template Walk<Hpolytope, RNGType> WalkType;
                    typedef RandomPointGenerator<WalkType> PointGenerator;
                    
                    if (nburns > 0) {
                        PointGenerator::apply(P, p, nburns, walk_length, rand_points,
                                            push_back_policy, rng);
                        rand_points.clear();
                    }
                    PointGenerator::apply(P, p, numpoints, walk_length, rand_points,
                                        push_back_policy, rng);
                }
                break;
            case BALL_WALK:
                {
                    typedef typename BallWalk::template Walk<Hpolytope, RNGType> WalkType;
                    typedef RandomPointGenerator<WalkType> PointGenerator;
                    
                    if (nburns > 0) {
                        PointGenerator::apply(P, p, nburns, walk_length, rand_points,
                                            push_back_policy, rng);
                        rand_points.clear();
                    }
                    PointGenerator::apply(P, p, numpoints, walk_length, rand_points,
                                        push_back_policy, rng);
                }
                break;
            default:
                error("sample_points: Unknown walk type");
                return octave_value_list();
            }
        }
        catch (const std::exception& e)
        {
            error("sample_points: Sampling failed: %s", e.what());
            return octave_value_list();
        }
        
        // Convert list to matrix (d x n, column-wise)
        result_matrix.resize(n, numpoints);
        unsigned int col = 0;
        for (typename std::list<Point>::iterator it = rand_points.begin(); 
             it != rand_points.end(); ++it, ++col)
        {
            for (int row = 0; row < n; ++row)
            {
                result_matrix(row, col) = (*it)[row];
            }
        }
        
        if (verbose)
        {
            octave_stdout << "[Volesti C++] Sampling complete!" << std::endl;
        }
    }
    else  // V-polytope
    {
        Matrix octave_V;
        int m, n;
        
        try
        {
            octave_V = args(0).matrix_value();
            m = octave_V.rows();  // number of vertices
            n = octave_V.cols();  // dimension
        }
        catch (const std::exception& e)
        {
            error("sample_points: Invalid input. Expected vertex matrix V. Error: %s", e.what());
            return octave_value_list();
        }
        
        if (m < n + 1)
        {
            error("sample_points: Need at least %d vertices for %dD V-polytope (got %d)", 
                  n+1, n, m);
            return octave_value_list();
        }
        
        // Auto-select walk type for V-polytope
        if (walk_type == -1)
        {
            walk_type = RDHR;  // Default for V-polytope
        }
        
        // Map vertex matrix to Eigen (zero-copy)
        Eigen::Map<MT> V_eigen(octave_V.fortran_vec(), m, n);
        
        // Create b vector (all ones)
        VT b_vec = VT::Ones(m);
        
        Vpolytope P(n, V_eigen, b_vec);
        
        if (verbose)
        {
            octave_stdout << "[Volesti C++] Sampling " << numpoints << " points from " 
                          << n << "D V-polytope (" << m << " vertices)..." << std::endl;
            const char* walk_names[] = {"CDHR", "RDHR", "BallWalk"};
            octave_stdout << "[Volesti C++] Walk: " << walk_names[walk_type]
                          << ", walk_length=" << walk_length 
                          << ", nburns=" << nburns << std::endl;
        }
        
        // Compute inner ball for starting point
        std::pair<Point, NT> inner_ball;
        try
        {
            inner_ball = P.ComputeInnerBall();
            if (inner_ball.second < 0.0)
            {
                error("sample_points: Unable to compute a feasible starting point");
                return octave_value_list();
            }
        }
        catch (const std::exception& e)
        {
            error("sample_points: Failed to compute inner ball: %s", e.what());
            return octave_value_list();
        }
        
        Point starting_point = inner_ball.first;
        
        // Create RNG
        RNGType rng(n);
        
        // Sample points using inline RandomPointGenerator pattern
        std::list<Point> rand_points;
        PushBackWalkPolicy push_back_policy;
        Point p = starting_point;
        
        try
        {
            switch (walk_type)
            {
            case CDHR:
                {
                    typedef typename CDHRWalk::template Walk<Vpolytope, RNGType> WalkType;
                    typedef RandomPointGenerator<WalkType> PointGenerator;
                    
                    if (nburns > 0) {
                        PointGenerator::apply(P, p, nburns, walk_length, rand_points,
                                            push_back_policy, rng);
                        rand_points.clear();
                    }
                    PointGenerator::apply(P, p, numpoints, walk_length, rand_points,
                                        push_back_policy, rng);
                }
                break;
            case RDHR:
                {
                    typedef typename RDHRWalk::template Walk<Vpolytope, RNGType> WalkType;
                    typedef RandomPointGenerator<WalkType> PointGenerator;
                    
                    if (nburns > 0) {
                        PointGenerator::apply(P, p, nburns, walk_length, rand_points,
                                            push_back_policy, rng);
                        rand_points.clear();
                    }
                    PointGenerator::apply(P, p, numpoints, walk_length, rand_points,
                                        push_back_policy, rng);
                }
                break;
            case BALL_WALK:
                {
                    typedef typename BallWalk::template Walk<Vpolytope, RNGType> WalkType;
                    typedef RandomPointGenerator<WalkType> PointGenerator;
                    
                    if (nburns > 0) {
                        PointGenerator::apply(P, p, nburns, walk_length, rand_points,
                                            push_back_policy, rng);
                        rand_points.clear();
                    }
                    PointGenerator::apply(P, p, numpoints, walk_length, rand_points,
                                        push_back_policy, rng);
                }
                break;
            default:
                error("sample_points: Unknown walk type");
                return octave_value_list();
            }
        }
        catch (const std::exception& e)
        {
            error("sample_points: Sampling failed: %s", e.what());
            return octave_value_list();
        }
        
        // Convert list to matrix (d x n, column-wise)
        result_matrix.resize(n, numpoints);
        unsigned int col = 0;
        for (typename std::list<Point>::iterator it = rand_points.begin(); 
             it != rand_points.end(); ++it, ++col)
        {
            for (int row = 0; row < n; ++row)
            {
                result_matrix(row, col) = (*it)[row];
            }
        }
        
        if (verbose)
        {
            octave_stdout << "[Volesti C++] Sampling complete!" << std::endl;
        }
    }

    return octave_value(result_matrix);
}
