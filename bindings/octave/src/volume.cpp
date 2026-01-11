/*
 * volume.cpp
 * 
 * GNU Octave interface for Volesti library
 * Volume computation for H-polytopes and V-polytopes
 * 
 * This wrapper demonstrates zero-copy architecture using Eigen::Map
 * to directly map Octave's memory into Volesti's Eigen structures.
 */

#include <octave/oct.h>
#include <octave/parse.h>
#include <Eigen/Eigen>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "volume/volume_sequence_of_balls.hpp"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;
typedef VPolytope<Point> Vpolytope;

DEFUN_DLD(compute_volume, args, nargout,
          "volume = compute_volume(A, b [, epsilon, walk_length, verbose])\n\
          volume = compute_volume(V [, epsilon, walk_length, verbose])\n\
\n\
Compute the volume of a polytope.\n\
\n\
H-Polytope (Ax <= b):\n\
  A : m x n matrix of constraint coefficients\n\
  b : m x 1 vector of constraint bounds\n\
\n\
V-Polytope (convex hull of vertices):\n\
  V : m x n matrix where each row is a vertex\n\
\n\
Optional parameters (both types):\n\
  epsilon : error tolerance (default: 1.0)\n\
  walk_length : random walk length (default: 1)\n\
  verbose : show progress messages (default: true)\n")
{
    if (args.length() < 1 || args.length() > 5)
    {
        error("compute_volume: 1 to 5 arguments required");
        return octave_value_list();
    }

    // Determine polytope type based on second argument
    // H-polytope: (A, b, ...) where b is column vector  
    // V-polytope: (V) or (V, epsilon, ...) where epsilon is scalar
    bool is_hpolytope = false;
    bool is_vpolytope = false;
    
    if (args.length() == 1)
    {
        is_vpolytope = true;
    }
    else if (args.length() >= 2)
    {
        // Check if second argument is scalar (V-poly) or vector (H-poly)
        if (args(1).is_scalar_type())
        {
            is_vpolytope = true;
        }
        else
        {
            is_hpolytope = true;
        }
    }

    // Parse optional parameters (shifted index for V-polytope)
    int param_offset = is_hpolytope ? 2 : 1;
    NT epsilon = 1.0;
    unsigned int walk_length = 1;
    bool verbose = true;
    
    if (args.length() >= param_offset + 1)
    {
        epsilon = args(param_offset).scalar_value();
        if (epsilon <= 0)
        {
            error("compute_volume: epsilon must be positive");
            return octave_value_list();
        }
    }
    
    if (args.length() >= param_offset + 2)
    {
        double walk_length_dbl = args(param_offset + 1).scalar_value();
        if (walk_length_dbl < 1 || walk_length_dbl > 1e6)
        {
            error("compute_volume: walk_length must be between 1 and 1e6");
            return octave_value_list();
        }
        walk_length = static_cast<unsigned int>(walk_length_dbl);
    }
    
    if (args.length() >= param_offset + 3)
    {
        verbose = args(param_offset + 2).bool_value();
    }

    NT volume = 0.0;
    
    if (is_hpolytope)
    {
        // H-Polytope path
        Matrix octave_A = args(0).matrix_value();
        ColumnVector octave_b = args(1).column_vector_value();
        
        int m = octave_A.rows();
        int n = octave_A.cols();
        
        if (octave_b.numel() != m)
        {
            error("compute_volume: Dimensions of A and b do not match");
            return octave_value_list();
        }
        
        // Zero-copy: Map Octave data to Eigen
        Eigen::Map<MT> A_eigen(octave_A.fortran_vec(), m, n);
        Eigen::Map<VT> b_eigen(octave_b.fortran_vec(), m);
        
        Hpolytope P(n, A_eigen, b_eigen);
        
        if (verbose)
        {
            octave_stdout << "[Volesti C++] Computing H-polytope volume (" << n << "D, "
                          << m << " constraints)..." << std::endl;
            octave_stdout << "[Volesti C++] epsilon=" << epsilon 
                          << ", walk_length=" << walk_length << std::endl;
        }
        
        try
        {
            volume = volume_sequence_of_balls(P, epsilon, walk_length);
        }
        catch (const std::exception& e)
        {
            error("compute_volume: H-polytope volume computation failed: %s", e.what());
            return octave_value_list();
        }
    }
    else  // V-polytope
    {
        Matrix octave_V = args(0).matrix_value();
        
        int m = octave_V.rows();  // number of vertices
        int n = octave_V.cols();  // dimension
        
        if (m < n + 1)
        {
            error("compute_volume: Need at least %d vertices for %dD V-polytope (got %d)", 
                  n+1, n, m);
            return octave_value_list();
        }
        
        // Zero-copy: Map vertex matrix to Eigen
        Eigen::Map<MT> V_eigen(octave_V.fortran_vec(), m, n);
        
        // Create b vector (all ones for standard V-polytope)
        VT b_vec = VT::Ones(m);
        
        Vpolytope P(n, V_eigen, b_vec);
        
        if (verbose)
        {
            octave_stdout << "[Volesti C++] Computing V-polytope volume (" << n << "D, "
                          << m << " vertices)..." << std::endl;
            octave_stdout << "[Volesti C++] epsilon=" << epsilon 
                          << ", walk_length=" << walk_length << std::endl;
        }
        
        try
        {
            volume = volume_sequence_of_balls(P, epsilon, walk_length);
        }
        catch (const std::exception& e)
        {
            error("compute_volume: V-polytope volume computation failed: %s", e.what());
            return octave_value_list();
        }
    }
    
    if (verbose)
    {
        octave_stdout << "[Volesti C++] Computation complete!" << std::endl;
    }

    return octave_value(volume);
}

