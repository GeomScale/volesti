function [ts, iP, iQ] = isAdjacent(P,Q,fP,fQ)
%
%  ISADJACENT: Test if a polyhedron shares a facet with another polyhedron. 
%  =========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = P.isAdjacent(Q)
%      ts = isAdjacent(P,Q)
%      [ts, iP, iQ] = isAdjacent(P,Q,fP,fQ)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the polyhedron P  has a facet to facet property with the
%  polyhedron Q. Both polyhedrons must be in H-representation. If they are not, the
%  irredundant H-representation will be computed. Basically, the function tests if
%  polyhedra P  and S  are adjacent by solving LP problem consecutively for each
%  facet. The polyhedra are declared as adjacent if their intersection is of
%  dimension d-1  and if the facet of polyhedron P  is also a facet for the
%  polyhedron S. If you want to test just specific facets, you can provide them in
%  fP and fQ arguments.
%  
%  INPUT
%  -----
%     
%        
%          P  Polyhedron in H-representation           
%             Class: Polyhedron                        
%          Q  Polyhedron in H-representation           
%             Class: Polyhedron                        
%          fP Index of a facet to test from polyhedron 
%             P.                                       
%             Class: double                            
%          fQ Index of a facet to test from polyhedron 
%             Q.                                       
%             Class: double                            
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts Logical statement if the polyhedron P    
%             is in a face to face property with Q.    
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%          iP Index of a facet from polyhedron P  that 
%             is common with polyhedron Q.             
%             Class: double                            
%          iQ Index of a facet from polyhedron Q  that 
%             is common with polyhedron P.             
%             Class: double                            
%               
%  
%  
%  SEE ALSO
%  --------
%     isNeighbor
%  

%  AUTHOR(s)
%  ---------
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

narginchk(2, 4);

if nargin < 3
    [ts, iP, iQ] = isNeighbor(P,Q);
elseif nargin < 4
    [ts, iP, iQ] = isNeighbor(P,Q,fP);
else
    [ts, iP, iQ] = isNeighbor(P,Q,fP,fQ);
end

% iP and iQ can have more elements which means that P or Q have
% equalities written as double-sided inequalities. In this case we take
% only the first index and that (although matlab seems to have no problem
% with checking multiple facets).
% see -> test_polyhedron_isadjacent_09_pass.m

if ts
    % check if facet of P and facet of Q are equal
    if ~isempty(iP)
        f1 = P.getFacet(iP(1));
    else 
        f1 = P;
    end
    if ~isempty(iQ)
        f2 = Q.getFacet(iQ(1));
    else
        f2 = Q;
    end
    if f1~=f2
        ts = false;
    end
end

end
 
