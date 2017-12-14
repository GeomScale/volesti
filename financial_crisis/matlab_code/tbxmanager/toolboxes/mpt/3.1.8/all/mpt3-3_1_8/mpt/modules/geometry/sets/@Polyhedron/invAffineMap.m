function Pnew = invAffineMap(P, T, t)
%
%  INVAFFINEMAP: Compute the inverse affine map of the Polyhedron. 
%  ================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = P.invAffineMap(T)
%      Q = P.invAffineMap(T, t)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes an inverse affine map Q of polyhedron Pbased on the transformation
%  matrix T  and vector t. The polyhedron Q is given by 
%                                    n                        
%                        Q = { x in R   |  Tx+t in P }     (1)
%                                                             
%     The matrix T must be a square real matrix. The vector t, if omitted, defaults
%  to a zero vector of corresponding dimension.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format.                
%            Class: Polyhedron                        
%          T Transformation matrix.                   
%            Class: double                            
%          t Transformation vector.                   
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q Polyhedron representing the affine map   
%            in H-representation.                     
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     affineHull
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
 
 
validate_realmatrix(T);
if nargin<3
	t = zeros(P(1).Dim, 1);
else
    validate_realvector(t);
end

% deal with arrays
if numel(P)>1
	Pnew(size(P)) = Polyhedron;
    for i=1:numel(P)
        Pnew(i) = P(i).invAffineMap(T, t);
    end
    return
end


% TODO: support non-square mappings
if size(T, 1)~=size(T, 2)
	error('Only square mappings supported.');
end

% TODO: deal with the V-representation directly
if ~P.hasHRep
	% we require the H-representation
	P.minHRep();
end

if isempty(P.He_int)
	% faster call if we have no equalities
	Pnew = Polyhedron(P.A*T, P.b-P.A*t);
else
	Pnew = Polyhedron('A', P.A*T, 'b', P.b-P.A*t, ...
		'Ae', P.Ae*T, 'be', P.be - P.Ae*t);
end

end
