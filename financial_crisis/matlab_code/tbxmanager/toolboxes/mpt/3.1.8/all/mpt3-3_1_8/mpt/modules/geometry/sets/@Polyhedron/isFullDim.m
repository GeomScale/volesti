function tf = isFullDim(P)
%
%  ISFULLDIM: Test if a polyhedron has a non-empty interior. 
%  ==========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      tf = P.isFullDim
%      tf = isFullDim(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the polyhedron P  has a non-empty interior and false
%  otherwise. A polyhedron has a non-empty interior if and only if its dimension is
%  the same as its dimension of representation. i.e. rank affhull P = size(P.A,2),
%  or equivalently if there exists a non-empty ball of dimension size(P.A,2) that
%  is contained within it.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          tf true if the polyhedron P  has a          
%             non-empty interior and false otherwise.  
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
%     isEmptySet,  isBounded
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
 
 
tf = false(size(P));
for i=1:length(P)
	% use a stored information if possible
	fulldim = P(i).Internal.FullDim;
	if ~isempty(fulldim)
		tf(i) = fulldim;
	else
		% if the polyhedron is empty -> not full dimensional
		if P(i).isEmptySet
			P(i).Internal.FullDim = false;
		else
			% compute interior point only if P(i).Empty property has not been set
			sol = P(i).interiorPoint;
			P(i).Internal.FullDim = sol.isStrict && ~isempty(sol.x);
		end
		tf(i) = P(i).Internal.FullDim;
	end
end
end
