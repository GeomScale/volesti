// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
//
// Octave interface for volume computation
//
// Licensed under GNU LGPL.3, see LICENSE file

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>
#include <octave/ov.h>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "preprocess/max_inscribed_ball.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef HPolytope<Point> HPOLYTOPE;

DEFUN_DLD (volesti_volume, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{vol} =} volesti_volume (@var{A}, @var{b})\n\
@deftypefnx {Loadable Function} {@var{vol} =} volesti_volume (@var{A}, @var{b}, @var{method})\n\
@deftypefnx {Loadable Function} {@var{vol} =} volesti_volume (@var{A}, @var{b}, @var{method}, @var{error})\n\
\n\
Compute the volume of an H-polytope defined by @var{A}*x <= @var{b}.\n\
\n\
@var{A} is a matrix of size m x d (m constraints, d dimensions).\n\
@var{b} is a column vector of size m.\n\
@var{method} is a string specifying the method:\n\
  - \"sequence_of_balls\" (default)\n\
  - \"cooling_gaussians\"\n\
@var{error} is the relative error tolerance (default: 0.1).\n\
\n\
Returns the estimated volume.\n\
@end deftypefn")
{
  int nargin = args.length ();
  
  if (nargin < 2 || nargin > 4)
    {
      error ("volesti_volume: wrong number of arguments");
      return octave_value_list ();
    }
  
  if (! args(0).is_matrix_type () || ! args(1).is_matrix_type ())
    {
      error ("volesti_volume: A and b must be matrices");
      return octave_value_list ();
    }
  
  // Get input matrices
  Matrix A_mat = args(0).matrix_value ();
  ColumnVector b_vec = args(1).column_vector_value ();
  
  if (A_mat.rows () != b_vec.length ())
    {
      error ("volesti_volume: number of rows in A must match length of b");
      return octave_value_list ();
    }
  
  // Convert to Eigen matrices
  unsigned int m = A_mat.rows ();
  unsigned int d = A_mat.cols ();
  
  Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> A(m, d);
  Eigen::Matrix<NT, Eigen::Dynamic, 1> b(m);
  
  for (unsigned int i = 0; i < m; i++)
    {
      b(i) = b_vec(i);
      for (unsigned int j = 0; j < d; j++)
        {
          A(i, j) = A_mat(i, j);
        }
    }
  
  // Create polytope
  HPOLYTOPE P(d, A, b);
  
  // Get method and error parameters
  std::string method = "sequence_of_balls";
  NT error = 0.1;
  
  if (nargin >= 3)
    {
      if (args(2).is_string ())
        method = args(2).string_value ();
    }
  
  if (nargin >= 4)
    {
      error = args(3).double_value ();
    }
  
  // Compute volume
  NT volume = 0.0;
  RNGType rng(d);
  
  try
    {
      if (method == "cooling_gaussians")
        {
          unsigned int walk_len = 10 + d / 10;
          volume = volume_cooling_gaussians<GaussianBallWalk, RNGType> (P, rng, error, walk_len);
        }
      else // sequence_of_balls (default)
        {
          unsigned int walk_len = 10 + d / 10;
          volume = volume_sequence_of_balls<CDHRWalk, RNGType> (P, rng, error, walk_len);
        }
    }
  catch (const std::exception& e)
    {
      error ("volesti_volume: error computing volume: %s", e.what ());
      return octave_value_list ();
    }
  
  return octave_value (volume);
}

