function [fval, feasible] = feval(obj, x, function_name)
%
%  FEVAL: Evaluates a function defined over a convex set or an array thereof. 
%  ===========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      [fval, feasible] = Set.feval(x)
%      [fval, feasible] = Set.feval(x,function_name)
%      [fval, feasible] = feval(Set,x,function_name)
%    
%  
%  DESCRIPTION
%  -----------
%     Evaluates function for given value of the point x over the convex Set
%  characterized by the the string function_name. The dimension of the vector x
%  must be the same as the dimension of the Set. If the Set is an array of convex
%  sets, the function values are returned in an array. For a point that lies
%  outside of the Set, the output is NaN.
%  
%  INPUT
%  -----
%     
%        
%          Set           Convex set or an array thereof, i.e. any 
%                        object derived from the ConvexSet class, 
%                        e.g. Polyhedron, YSet, ...               
%                        Class: ConvexSet                         
%          x             A point at which the function should be  
%                        evaluated. The point must be given as    
%                        column and must be in the same dimension 
%                        as the set.                              
%                        Class: double                            
%          function_name String name of the function to evaluate. 
%                        It must refer to a single function. If   
%                        omitted, S.feval(x) only works if the    
%                        set has a single function.               
%                        Class: char                              
%                          
%  
%  
%  OUTPUT
%  ------
%     
%        
%          fval     An (n X m)  matrix of function values at 
%                   x, where m  is the number of sets. If x  
%                   is not contained in the j-th set, then   
%                   the j-th column of fval is a vector of   
%                   NaN.                                     
%                   Class: double                            
%          feasible An  (1 X m)  vector of logicals.         
%                   feasible(j)=true if the j-th element of  
%                   the array contains x.                    
%                   Class: double                            
%                     
%  
%  
%  SEE ALSO
%  --------
%     fplot,  Function
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
%  
%  

%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%    This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
narginchk(2, 3);
if nargin==2
	function_name = '';
end
n_obj = numel(obj);
if n_obj==0
	% exit quickly if we have an empty set (in the matlab sense of empty)
	fval = [];
	feasible = false;
	return
end

%% validate arguments
[function_name, msg] = obj.validateFunctionName(function_name);
error(msg); % the error is only thrown if msg is not empty

%% evaluate
validate_realvector(x);
if n_obj==1
	% faster implementation for single sets
	feasible = obj.contains(x);
	fval = obj.Functions(function_name).feval(x);
	if ~feasible
		fval = NaN(size(fval));
	end
else
	% arrays of sets
	feasible = false(1, n_obj);
	fval = [];
	for i = 1:n_obj
		feasible(i) = obj(i).contains(x);
		f = obj(i).Functions(function_name).feval(x);
		if ~feasible(i)
			f = NaN(size(f));
		end
		fval = [fval f];
	end
end

end
