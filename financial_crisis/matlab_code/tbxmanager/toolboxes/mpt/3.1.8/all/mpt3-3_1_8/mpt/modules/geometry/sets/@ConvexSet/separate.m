function sep = separate(obj, x)
%
%  SEPARATE: Computes separating hyperplane between the set and given point. 
%  ==========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      h = separate(S,x)
%      h = S.separate(x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute a separating hyperplane between the convex Set and the given point x.
%  If Set is an array, the hyperplanes are returned as a cell array.
%  
%  INPUT
%  -----
%     
%        
%          S Any set derived from ConvexSet class,    
%            e.g. YSet or Polyhedron.                 
%            Class: ConvexSet                         
%          x The point given as vector in the same    
%            dimension as ConvexSet.                  
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          h Separating hyperplane defined as {  x    
%              h(                                     
%             x                                       
%              -1 )=0 }.                              
%            Class: double                            
%              
%  
%  
%  SEE ALSO
%  --------
%     distance,  support
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
 
 
global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);
error(obj.rejectArray());

validate_realvector(x);

if numel(x)~=obj.Dim
    error('The point "x" must be a vector with the length of %i.', obj.Dim);
end

% for empty objects, quickly return
if obj.isEmptySet
    sep = [];
    return;
end

x = x(:);
ret = obj.project(x);
switch ret.exitflag
    case MPTOPTIONS.INFEASIBLE,
        sep = [];
        return
    case MPTOPTIONS.UNBOUNDED,
        error('Solver returned unbounded. This should not happen');
end
if ret.dist < MPTOPTIONS.rel_tol,
    sep = [];
    return;
end

y = ret.x; % Projected point

sep = y-x; % Hyperplane normal points from x to y
sep(end+1) = sep'*(y+x)/2; % Hyperplane intersects point halfway between x and y
sep = sep(:)';

end
