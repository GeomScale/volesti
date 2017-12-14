function status = le(U1, U2)
%
%  LE: Test if a union of polyhedra is contained inside another union. 
%  ====================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U1 <= U2
%    
%  
%  DESCRIPTION
%  -----------
%     Check if the union of polyhedra U1 is a non-strict subset of the union U2.
%  The result it the logical statement if U1 <= U2 and false otherwise.
%  
%  INPUT
%  -----
%     
%        
%          U1 Union of polyhedra in the same           
%             dimension.                               
%             Class: PolyUnion                         
%          U2 Union of polyhedra in the same           
%             dimension.                               
%             Class: PolyUnion                         
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          tf True if U1 <= U2 and false otherwise.    
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
%     ge
%  

%  AUTHOR(s)
%  ---------
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
 
 
global MPTOPTIONS

if isa(U1, 'Polyhedron')
	U1 = PolyUnion(U1);
end
if isa(U2, 'Polyhedron')
	U2 = PolyUnion(U2);
end
if ~isa(U1,'PolyUnion') || ~isa(U2,'PolyUnion')
    error('All inputs must be PolyUnion objects.');
end
% use forEach if you have multiple unions
error(U1.rejectArray());
error(U2.rejectArray());

if (U1.Num==0 || all(isEmptySet(U1.Set)))
	% empty set is contained in any set
	status = true;
	return
end
if (U2.Num==0 || all(isEmptySet(U2.Set)));
	% non-empty set U1 cannot contain an empty set
	status = false;
	return
end

% heuristics: check containement of outer approximations
B1 = U1.outerApprox();
B2 = U2.outerApprox();
bbox_tol = 10*MPTOPTIONS.rel_tol;
if any(B1.Internal.lb + bbox_tol < B2.Internal.lb) || ...
        any(B1.Internal.ub - bbox_tol > B2.Internal.ub),
    % outer approximation of B1 is not contained in the outer
    % approximation of B2, hence B1 cannot be contained in B2
    status = false;
    return
end

if U2.Num==1
	% simpler case, test whether each elements of U1 is contained in
	% U2.Set(1)
	for i = 1:U1.Num
		if ~(U1.Set(i) <= U2.Set(1))
			status = false;
			return
		end
	end
	status = true;

else
	% complicated case, compute set difference using the "noconstruction"
	% mode
	status = all(isEmptySet(mldivide(U1.Set, U2.Set, true)));
	
end

end
