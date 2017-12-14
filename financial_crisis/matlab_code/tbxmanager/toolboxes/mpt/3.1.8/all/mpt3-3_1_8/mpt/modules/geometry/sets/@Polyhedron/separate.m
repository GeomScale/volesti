function sep = separate(P, S)
%
%  SEPARATE: Separate a point/polyhedron from another polyhedron. 
%  ===============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      h = P.separate(x)
%      h = separate(P, x)
%      h = P.separate(S)
%      h = separate(P, S)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes a separating hyperplane h  between P  and x/ S  : 
%                       T (    )                     T (    )            
%           h = {y  |  h  ( y  )<= 0 forall y in P, h  ( x  )= 0 } >= 0 }
%                         ( -1 )                       ( -1 )            
%     or 
%                T (    )                     T (    )                          
%    h = {y  |  h  ( y  )<= 0 forall y in P, h  ( z  )= 0 } >= 0 forall z in S }
%                  ( -1 )                       ( -1 )                          
%  
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          S Polyhedron in any format                 
%            Class: Polyhedron                        
%          x Column vector of length P.Dim.           
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          h Separating hyperplane {x  |  h^T (       
%             x                                       
%             -1 ) = 0 }, or [] if none exists.       
%            Class: double                            
%              
%  
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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

if isnumeric(S)
    sep = P.separate@ConvexSet(S);
    return;
end

if ~isa(S, 'Polyhedron'),
    error('S must be a polyhedron.');
end
if numel(S)>1
    error('Only single polyhedron S is allowed.');
end

% deal with arrays
error(P.rejectArray());

if S.Dim ~= P.Dim,
    error('S must be in the dimension %i, the same as P.', P.Dim);
end

% Compute the closest points between P and S
ret = distance(P, S);

if ret.exitflag~=MPTOPTIONS.OK
    % infeasible
    sep = [];
    return;
end

% if the sets intersect
if ret.dist < MPTOPTIONS.rel_tol,
    sep = [];
    return;
end


% Compute a separating hyperplane
sep = ret.y-ret.x; % Hyperplane normal points from x to y
sep(end+1) = sep'*(ret.y+ret.x)/2; % Hyperplane intersects point halfway between x and y
sep = sep(:)';


end
