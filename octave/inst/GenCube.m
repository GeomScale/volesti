## -*- texinfo -*-
## @deftypefn {Function File} {@var{P} =} GenCube (@var{d})
## @deftypefnx {Function File} {@var{P} =} GenCube (@var{d}, @var{representation})
## @deftypefnx {Function File} {@var{P} =} GenCube (@var{d}, @var{representation}, @var{scale})
##
## Generate a hypercube polytope of dimension @var{d}.
##
## @var{d} is the dimension of the cube.
## @var{representation} is a string:
##   - "H" (default) - H-representation (Ax <= b)
##   - "V" - V-representation (not yet supported in Octave interface)
##
## @var{scale} is the scaling factor (default: 1.0).
##
## Returns a polytope structure with fields A and b.
##
## @seealso{volume, sample_points}
## @end deftypefn

function P = GenCube (d, representation = "H", scale = 1.0)
  if (nargin < 1)
    print_usage ();
  endif
  
  if (! isscalar (d) || d <= 0)
    error ("GenCube: d must be a positive integer");
  endif
  
  if (! ischar (representation))
    error ("GenCube: representation must be a string");
  endif
  
  if (strcmpi (representation, "V"))
    error ("GenCube: V-representation not yet supported in Octave interface");
  endif
  
  % Generate H-representation of cube: -scale <= x_i <= scale
  % This gives: x_i <= scale and -x_i <= scale
  % Which is: [I; -I] * x <= [scale; scale]
  A = [eye(d); -eye(d)];
  b = scale * ones(2*d, 1);
  
  P = struct ("A", A, "b", b, "dimension", d);
endfunction

