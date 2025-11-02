## -*- texinfo -*-
## @deftypefn {Function File} {@var{vol} =} volume (@var{P})
## @deftypefnx {Function File} {@var{vol} =} volume (@var{A}, @var{b})
## @deftypefnx {Function File} {@var{vol} =} volume (@var{A}, @var{b}, @var{method})
## @deftypefnx {Function File} {@var{vol} =} volume (@var{A}, @var{b}, @var{method}, @var{error})
##
## Compute the volume of a convex polytope.
##
## If called with a polytope structure @var{P}, compute its volume.
## If called with @var{A} and @var{b}, compute the volume of the H-polytope
## defined by @var{A}*x <= @var{b}.
##
## @var{method} is a string specifying the volume computation method:
##   - "sequence_of_balls" (default)
##   - "cooling_gaussians"
##
## @var{error} is the relative error tolerance (default: 0.1).
##
## @seealso{volesti_sample, GenCube}
## @end deftypefn

function vol = volume (varargin)
  if (nargin == 0)
    print_usage ();
  endif
  
  % Handle polytope structure
  if (nargin == 1 && isstruct (varargin{1}))
    P = varargin{1};
    if (! isfield (P, "A") || ! isfield (P, "b"))
      error ("volume: polytope structure must have fields A and b");
    endif
    A = P.A;
    b = P.b;
    method = "sequence_of_balls";
    error = 0.1;
  elseif (nargin >= 2)
    A = varargin{1};
    b = varargin{2};
    method = "sequence_of_balls";
    error = 0.1;
    
    if (nargin >= 3)
      method = varargin{3};
    endif
    if (nargin >= 4)
      error = varargin{4};
    endif
  else
    print_usage ();
  endif
  
  % Validate inputs
  if (! ismatrix (A) || ! ismatrix (b))
    error ("volume: A and b must be matrices");
  endif
  
  if (rows (A) != rows (b))
    error ("volume: number of rows in A must match length of b");
  endif
  
  % Call MEX function
  vol = volesti_volume (A, b, method, error);
endfunction

