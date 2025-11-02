// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
//
// Octave interface for sampling from polytopes
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
#include "sampling/sampling.hpp"
#include "sampling/random_point_generators.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "preprocess/max_inscribed_ball.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef HPolytope<Point> HPOLYTOPE;
typedef std::list<Point> PointList;

DEFUN_DLD (volesti_sample, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{points} =} volesti_sample (@var{A}, @var{b}, @var{n})\n\
@deftypefnx {Loadable Function} {@var{points} =} volesti_sample (@var{A}, @var{b}, @var{n}, @var{method})\n\
@deftypefnx {Loadable Function} {@var{points} =} volesti_sample (@var{A}, @var{b}, @var{n}, @var{method}, @var{walk_len})\n\
@deftypefnx {Loadable Function} {@var{points} =} volesti_sample (@var{A}, @var{b}, @var{n}, @var{method}, @var{walk_len}, @var{nburns})\n\
\n\
Sample uniform points from an H-polytope defined by @var{A}*x <= @var{b}.\n\
\n\
@var{A} is a matrix of size m x d (m constraints, d dimensions).\n\
@var{b} is a column vector of size m.\n\
@var{n} is the number of points to sample.\n\
@var{method} is a string specifying the random walk method:\n\
  - \"cdhr\" (default) - Coordinate Directions Hit-and-Run\n\
  - \"rdhr\" - Random Directions Hit-and-Run\n\
  - \"ball\" - Ball Walk\n\
  - \"billiard\" - Billiard Walk\n\
@var{walk_len} is the walk length (default: 10 + dimension/10).\n\
@var{nburns} is the number of burn-in steps (default: 0).\n\
\n\
Returns a d x n matrix of sampled points.\n\
@end deftypefn")
{
  int nargin = args.length ();
  
  if (nargin < 3 || nargin > 6)
    {
      error ("volesti_sample: wrong number of arguments");
      return octave_value_list ();
    }
  
  if (! args(0).is_matrix_type () || ! args(1).is_matrix_type ())
    {
      error ("volesti_sample: A and b must be matrices");
      return octave_value_list ();
    }
  
  if (! args(2).is_scalar_type () || args(2).double_value () <= 0)
    {
      error ("volesti_sample: n must be a positive integer");
      return octave_value_list ();
    }
  
  // Get input matrices
  Matrix A_mat = args(0).matrix_value ();
  ColumnVector b_vec = args(1).column_vector_value ();
  unsigned int n = static_cast<unsigned int> (args(2).int_value ());
  
  if (A_mat.rows () != b_vec.length ())
    {
      error ("volesti_sample: number of rows in A must match length of b");
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
  
  // Get starting point (Chebychev center or origin)
  auto InnerBall = P.ComputeInnerBall();
  Point starting_point;
  if (InnerBall.second > 0.0)
    {
      starting_point = InnerBall.first;
    }
  else
    {
      // Use origin as fallback
      starting_point = Point(d);
    }
  
  // Get method and parameters
  std::string method = "cdhr";
  unsigned int walk_len = 10 + d / 10;
  unsigned int nburns = 0;
  
  if (nargin >= 4)
    {
      if (args(3).is_string ())
        method = args(3).string_value ();
    }
  
  if (nargin >= 5)
    {
      walk_len = static_cast<unsigned int> (args(4).int_value ());
    }
  
  if (nargin >= 6)
    {
      nburns = static_cast<unsigned int> (args(5).int_value ());
    }
  
  // Sample points
  PointList randPoints;
  RNGType rng(d);
  
  try
    {
      if (method == "rdhr")
        {
          uniform_sampling<RDHRWalk> (randPoints, P, rng, walk_len, n, starting_point, nburns);
        }
      else if (method == "ball")
        {
          uniform_sampling<BallWalk> (randPoints, P, rng, walk_len, n, starting_point, nburns);
        }
      else if (method == "billiard")
        {
          uniform_sampling<BilliardWalk> (randPoints, P, rng, walk_len, n, starting_point, nburns);
        }
      else // cdhr (default)
        {
          uniform_sampling<CDHRWalk> (randPoints, P, rng, walk_len, n, starting_point, nburns);
        }
    }
  catch (const std::exception& e)
    {
      error ("volesti_sample: error sampling points: %s", e.what ());
      return octave_value_list ();
    }
  
  // Convert points to Octave matrix
  Matrix points(d, n);
  unsigned int col = 0;
  for (auto pit = randPoints.begin (); pit != randPoints.end (); ++pit, ++col)
    {
      if (col >= n) break;
      auto coords = pit->getCoefficients ();
      for (unsigned int row = 0; row < d; row++)
        {
          points(row, col) = coords(row);
        }
    }
  
  return octave_value (points);
}

