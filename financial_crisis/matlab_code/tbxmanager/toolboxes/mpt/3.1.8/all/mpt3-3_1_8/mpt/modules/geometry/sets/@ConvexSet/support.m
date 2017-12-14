function supp = support(obj, x)
%
%  SUPPORT: Compute the support of the set in the specified direction. 
%  ====================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      s = support(S,x)
%      s = S.support(x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the support of the set in the direction given by the point x. The
%  underlying optimization problem to be solved is as follows 
%                                       ymax   x'y            
%                                 s.t.    y in Set            
%     where x  is the point with the desired direction and y  is the point lying
%  inside the convex Set. The support is returned as the optimal value of the
%  objective function  x'y. The dimension of x must be the same as the Set. If an
%  error occurs during by solving the above optimization problem, the support is
%  empty.
%  
%  INPUT
%  -----
%     
%        
%          S Any set derived from ConvexSet class,    
%            e.g. YSet or Polyhedron.                 
%            Class: ConvexSet                         
%          x The point given as real vector in the    
%            same dimension as the ConvexSet.         
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          s The support is returned as the optimal   
%            value of the cost function.              
%            Class: double                            
%              
%  
%  
%  SEE ALSO
%  --------
%     separate,  distance
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

% x can be an array of points put in a matrix
validate_realmatrix(x);

% deal with arrays
if numel(obj)>1
	supp = obj.forEach(@(e) e.support(x));
	return
end

if ~isequal(size(x, 1), obj.Dim)
	error('Input argument "x" must have %d rows.', obj.Dim);
end

% compute the support for each point (points are assumed to be stored
% column-wise)
n_points = size(x, 2);
supp = Inf*ones(n_points, 1);
for i=1:n_points
    sol = obj.extreme(x(:,i));
    if ~isempty(sol.supp)
        supp(i) = sol.supp;
    end
end

end
