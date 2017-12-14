function ret = distance(obj, x)
%
%  DISTANCE: Computes the closest distance between the convex set and given point. 
%  ================================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ret = distance(Set,x)
%      ret = Set.distance(x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the closest distance between the convex Set and given point xThe
%  approach is based on solving the optimization problem 
%                                                             
%                                  ymin   ||x-y||             
%                                                2            
%                                 s.t.    y in Set            
%     where x  is the given point in the same dimension as the set and y  is the
%  point inside the Set. If the optimization terminated successfully, the output
%  contains the distance ret.d and the points ret.y, ret.x. Otherwise, the output
%  is empty. If the Set is an array of convex sets, the distance and the point y
%  are returned in a cell arrays.
%  
%  INPUT
%  -----
%     
%        
%          Set Any object derived from the ConvexSet    
%              class, e.g. Polyhedron, YSet, ...        
%              Class: ConvexSet                         
%          x   A point given as a real vector with the  
%              same dimension as the convex set.        
%              Class: double                            
%                
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ret          Structure that contains the information  
%                       about the computed distance and the      
%                       points x, y.                             
%                       Class: struct                            
%          ret.exitflag Termination status from the related      
%                       optimization problem.                    
%                       Class: double                            
%          ret.dist     Distance between the point x and the     
%                       convex Set.                              
%                       Class: double                            
%          ret.y        The point that is contained in the       
%                       convex Set and is closest to x.          
%                       Class: double                            
%          ret.x        The point x that was provided.           
%                       Class: double                            
%                         
%  
%  
%  SEE ALSO
%  --------
%     outerApprox,  support,  separate
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
narginchk(2, 2);

% if obj is an array, put the results inside an array
if numel(obj)>1
	% return an array of structures
	ret = obj.forEach(@(elem) elem.distance(x));
    return
end

validate_realvector(x);

if length(x)~=obj.Dim
    error('The vector must have a length of %i', obj.Dim);
end

sol = obj.project(x);
ret.exitflag = sol.exitflag;
ret.dist = sol.dist;
ret.x = x;
ret.y = sol.x;
end
