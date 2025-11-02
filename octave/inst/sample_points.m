## -*- texinfo -*-
## @deftypefn {Function File} {@var{points} =} sample_points (@var{P}, @var{n})
## @deftypefnx {Function File} {@var{points} =} sample_points (@var{A}, @var{b}, @var{n})
## @deftypefnx {Function File} {@var{points} =} sample_points (@var{A}, @var{b}, @var{n}, @var{method})
## @deftypefnx {Function File} {@var{points} =} sample_points (@var{A}, @var{b}, @var{n}, @var{method}, @var{walk_len})
## @deftypefnx {Function File} {@var{points} =} sample_points (@var{A}, @var{b}, @var{n}, @var{method}, @var{walk_len}, @var{nburns})
##
## Sample uniform points from a convex polytope.
##
## If called with a polytope structure @var{P} and number of points @var{n},
## sample @var{n} points from the polytope.
## If called with @var{A}, @var{b}, and @var{n}, sample from the H-polytope
## defined by @var{A}*x <= @var{b}.
##
## @var{method} is a string specifying the random walk method:
##   - "cdhr" (default) - Coordinate Directions Hit-and-Run
##   - "rdhr" - Random Directions Hit-and-Run
##   - "ball" - Ball Walk
##   - "billiard" - Billiard Walk
##
## @var{walk_len} is the walk length (default: 10 + dimension/10).
## @var{nburns} is the number of burn-in steps (default: 0).
##
## Returns a d x n matrix of sampled points, where d is the dimension.
##
## @seealso{volume, GenCube}
## @end deftypefn

function points = sample_points (varargin)
  if (nargin < 2)
    print_usage ();
  endif
  
  % Handle polytope structure
  if (nargin == 2 && isstruct (varargin{1}))
    P = varargin{1};
    if (! isfield (P, "A") || ! isfield (P, "b"))
      error ("sample_points: polytope structure must have fields A and b");
    endif
    A = P.A;
    b = P.b;
    n = varargin{2};
    method = "cdhr";
    walk_len = [];
    nburns = 0;
  elseif (nargin >= 3)
    A = varargin{1};
    b = varargin{2};
    n = varargin{3};
    method = "cdhr";
    walk_len = [];
    nburns = 0;
    
    if (nargin >= 4)
      method = varargin{4};
    endif
    if (nargin >= 5)
      walk_len = varargin{5};
    endif
    if (nargin >= 6)
      nburns = varargin{6};
    endif
  else
    print_usage ();
  endif
  
  % Validate inputs
  if (! ismatrix (A) || ! ismatrix (b))
    error ("sample_points: A and b must be matrices");
  endif
  
  if (rows (A) != rows (b))
    error ("sample_points: number of rows in A must match length of b");
  endif
  
  if (! isscalar (n) || n <= 0)
    error ("sample_points: n must be a positive integer");
  endif
  
  % Call MEX function with appropriate arguments
  if (isempty (walk_len))
    points = volesti_sample (A, b, n, method);
  else
    points = volesti_sample (A, b, n, method, walk_len, nburns);
  endif
endfunction

