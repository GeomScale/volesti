function sol = facetInteriorPoints(P)
%
%  FACETINTERIORPOINTS: Compute points that lie on each of the facet of the
%  ========================================================================
%  Polyhedron. 
%  ============
%  
%  
%  SYNTAX
%  ------
%     
%      x = P.facetInteriorPoints
%      x = facetInteriorPoints(P)
%    
%  
%  DESCRIPTION
%  -----------
%     For each of the facet of the polyhedron compute a point in the relative
%  interior of that facet. The output is a matrix formed by concatenating the
%  interior point row-wise. Note that the polyhedron must be in the minimal
%  H-representation in order to proceed.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format.                
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          x Matrix formed by interior points where   
%            the row of the matrix corresponds to     
%            facet of the polyhedron.                 
%            Class: double                            
%              
%  
%  
%  SEE ALSO
%  --------
%     chebyCenter,  interiorPoint
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

% use P.forEach() for arrays
error(P.rejectArray());

%% compute interior points via projection
% check if P is not empty
if P.isEmptySet
    sol = [];
    return;
end

% Can't compute points in the facets without the facets
P.minHRep();

% Compute interior point of the polyhedron
xp = P.interiorPoint;
% get the point
x = xp.x;

m = size(P.H,1);
sol = zeros(m, P.Dim);
for j=1:m
    n = P.H(j,1:end-1);
    
    % Project x onto each facet to see if it's interior
    alpha = (P.H(j,end) - n*x) / (n*n');
    y = x + alpha * n';
    
    % Test if y is in the interior of the facet
    s = P.H*[y;-1];
    s(j) = [];
    % test equalities as well
    t = P.He*[y;-1];
    if all(s < MPTOPTIONS.abs_tol) && norm(t,Inf)<MPTOPTIONS.rel_tol
        sol(j,:) = y';
    else
        % Compute an interior point on the facet by solving an LP
        s = P.interiorPoint(j);
        if isempty(s.x)
            error('Could not compute interior point in the %i''th facet.', j);
        end
        sol(j,:) = s.x';
    end
    
end

end
