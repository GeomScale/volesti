/*
 * volesti_volume.cpp
 * 
 * GNU Octave interface for Volesti library
 * Proof of Concept: Volume computation for H-polytopes
 * 
 * This wrapper demonstrates zero-copy architecture using Eigen::Map
 * to directly map Octave's memory into Volesti's Eigen structures.
 */

#include <octave/oct.h>
#include <octave/parse.h>
#include <Eigen/Eigen>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "volume/volume_sequence_of_balls.hpp"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;

DEFUN_DLD(compute_volume, args, nargout,
          "volume = compute_volume(A, b [, epsilon, walk_length])\n\
\n\
Compute the volume of an H-polytope defined by Ax <= b.\n\
\n\
Parameters:\n\
  A : m x n matrix of constraint coefficients\n\
  b : m x 1 vector of constraint bounds\n\
  epsilon : (optional) error tolerance for approximation (default: 1.0)\n\
            Smaller values give more accurate results but take longer\n\
  walk_length : (optional) random walk length (default: 1)\n\
                Larger values improve mixing but increase computation time\n\
\n\
Returns:\n\
  volume : Estimated volume of the polytope\n\
\n\
Examples:\n\
  % Basic usage (default parameters)\n\
  A = [1 0; -1 0; 0 1; 0 -1];\n\
  b = ones(4, 1);\n\
  vol = compute_volume(A, b)\n\
\n\
  % Higher accuracy (epsilon = 0.1)\n\
  vol = compute_volume(A, b, 0.1)\n\
\n\
  % Custom epsilon and walk_length\n\
  vol = compute_volume(A, b, 0.1, 10)\n")
{
    if (args.length() < 2 || args.length() > 4)
    {
        error("compute_volume: 2 to 4 arguments required (A, b [, epsilon, walk_length])");
        return octave_value_list();
    }

    
    Matrix octave_A = args(0).matrix_value();
    ColumnVector octave_b = args(1).column_vector_value();


    int m = octave_A.rows();    // number of constraints
    int n = octave_A.cols();    // dimension of space


    if (octave_b.numel() != m)
    {
        error("compute_volume: Dimensions of A and b do not match");
        return octave_value_list();
    }

    if (n < 1 || m < n + 1)
    {
        error("compute_volume: Invalid polytope dimensions");
        return octave_value_list();
    }

    // Parse optional parameters
    NT epsilon = 1.0;  // Default error tolerance
    unsigned int walk_length = 1;  // Default walk length
    
    if (args.length() >= 3)
    {
        epsilon = args(2).scalar_value();
        if (epsilon <= 0)
        {
            error("compute_volume: epsilon must be positive");
            return octave_value_list();
        }
    }
    
    if (args.length() >= 4)
    {
        walk_length = static_cast<unsigned int>(args(3).scalar_value());
        if (walk_length < 1)
        {
            error("compute_volume: walk_length must be >= 1");
            return octave_value_list();
        }
    }

    // ZERO-COPY MAGIC: Map Octave data directly to Eigen structures
    // This avoids expensive memory duplication for large matrices
    Eigen::Map<MT> A_eigen(octave_A.fortran_vec(), m, n);
    Eigen::Map<VT> b_eigen(octave_b.fortran_vec(), m);

    
    Hpolytope P(n, A_eigen, b_eigen);

    // WATERMARK: Prove C++ execution 
    octave_stdout << "[Volesti C++] Computing volume for " << n << "D polytope with " 
                  << m << " constraints..." << std::endl;
    octave_stdout << "[Volesti C++] Parameters: epsilon=" << epsilon 
                  << ", walk_length=" << walk_length << std::endl;
    octave_stdout << "[Volesti C++] Using stochastic approximation (volume_sequence_of_balls)" 
                  << std::endl;

    
    NT volume = 0.0;
    
    try
    {
        volume = volume_sequence_of_balls(P, epsilon, walk_length);
        octave_stdout << "[Volesti C++] Computation complete!" << std::endl;
    }
    catch (const std::exception& e)
    {
        error("compute_volume: Volume computation failed: %s", e.what());
        return octave_value_list();
    }

    return octave_value(volume);
}

