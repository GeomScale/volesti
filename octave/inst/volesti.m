## -*- texinfo -*-
## @deftypefn {Function File} {} volesti ()
## @deftypefnx {Function File} {} volesti (@var{version})
## 
## Display information about the VolEsti package.
##
## If called with no arguments, display version and copyright information.
## If called with the argument "version", return the version string.
##
## @seealso{}
## @end deftypefn

function retval = volesti (varargin)
  if (nargin > 0 && strcmpi (varargin{1}, "version"))
    retval = "1.0.0";
    return;
  endif
  
  printf ("\n");
  printf ("VolEsti - Volume Estimation and Sampling of Convex Bodies\n");
  printf ("Version: 1.0.0\n");
  printf ("\n");
  printf ("VolEsti is a C++ library for volume approximation and sampling\n");
  printf ("of convex bodies (e.g. polytopes).\n");
  printf ("\n");
  printf ("This package provides an Octave interface to the VolEsti library.\n");
  printf ("\n");
  printf ("Copyright (c) 2012-2024 Vissarion Fisikopoulos\n");
  printf ("Copyright (c) 2018-2024 Apostolos Chalkis\n");
  printf ("Copyright (c) 2020-2024 Elias Tsigaridas\n");
  printf ("\n");
  printf ("Licensed under GNU LGPL.3\n");
  printf ("\n");
endfunction

