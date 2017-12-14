function ts = eq(P, S)
%
%  EQ: Returns true if the set covered by polyhedra P  is the same as the set
%  ==========================================================================
%  covered by S  and false otherwise. 
%  ===================================
%  
%  
%  SYNTAX
%  ------
%     
%      tf = P.eq(S)
%      tf = P == S
%    
%  
%  DESCRIPTION
%  -----------
%     Returns true if P  equals S  and false otherwise by testing if both P
%  subseteq S  and S subseteq P.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format or array of     
%            Polyhedra in H-representation.           
%            Class: Polyhedron                        
%          S Polyhedron (or array of polyhedra) in    
%            the same dimension as P.                 
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          tf True if S == P and false otherwise.      
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%               
%  
%  
%  SEE ALSO
%  --------
%     neq,  contains,  le,  lt,  ge,  gt
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
 
 
validate_polyhedron(S);

% both polyhedra are empty arrays
if numel(P)==0 && numel(S)==0
    ts = true;
    return
end
% of the polyhedra is empty array
if numel(P)==0 || numel(S)==0
    ts = false;
    return
end
    
dimP = [P.Dim];
if any(diff(dimP))
    error('All polyhedra "P" must be of the same dimension.');
end
dimS = [S.Dim];
if any(diff(dimS))
    error('All polyhedra "S" must be of the same dimension.');
end

if dimP(1)~=dimS(1)
  error('Polyhedra must be of the same dimension.');
end

if numel(P)>1 || numel(S)>1
	% first compare outer approximations
    if all(P.isEmptySet())
        % if all sets in P are empty, PolyUnion(P) produces an empty object
        % for which the outer approximation is always an empty set in R^0.
        % Thus we need to maintain dimensions. (issue #110)
        B1 = P(1);
    else
        B1 = PolyUnion(P).outerApprox;
    end
    if all(S.isEmptySet())
        B2 = S(1);
    else
        B2 = PolyUnion(S).outerApprox;
    end
	if ~(B1==B2)
		% bounding boxes differ, sets cannot be equal
		ts = false;
	else
		no_construction = true;
		ts = all(isEmptySet(mldivide(S, P, no_construction))) && ...
			all(isEmptySet(mldivide(P, S, no_construction)));
		% Note the special syntax of mldivide with three input arguments,
		% which makes checking for P\S==0 and S\P==0 much faster.
	end
else
    ts = P.contains(S) && S.contains(P);    
end

end
