function lcp_test
%
%  LCP_TEST: Example script for LCP solver 
%  ========================================
%  
%  
%  SYNTAX
%  ------
%     
%      lcp_test
%    
%  
%  DESCRIPTION
%  -----------
%     Runs a simple example to demonstrate LCP solver in Matlab and Simulink.
%  
%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%  
%  
%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under the terms of the GNU
%  General Public License as published by the Free Software Foundation; either version 2.1 of the
%  License, or (at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%  General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this library; if not,
%  write to the  Free Software Foundation, Inc.,  59 Temple Place, Suite 330,  Boston, MA 02111-1307
%  USA
 
 
 
 
 
 
disp('testing simple example');
fprintf(' M = [1 -1 -1 -1;\n      -1 1 -1 -1;\n     1 1 2 0;\n     1 1 0 2];\n')
disp('q = [3, 5, -9, -5]');
fprintf('syntax: \n     [z, w, basis, exitflag, pivots] = lcp(M, q, options)\n');
[z,w,basis,exitflag,pivots] = lcp([1 -1 -1 -1;-1 1 -1 -1; 1 1 2 0;1 1 0 2],[3;5;-9;-5])

disp('opening simulink interface')
lcp_sim

end
