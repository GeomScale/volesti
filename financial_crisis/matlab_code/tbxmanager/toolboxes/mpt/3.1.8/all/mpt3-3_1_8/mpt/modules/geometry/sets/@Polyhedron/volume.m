function vol = volume(P)
%
%  VOLUME: Compute the volume of the polyhedron. 
%  ==============================================
%  
%  
%  SYNTAX
%  ------
%     
%      vol = P.volume
%      vol = volume(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes the volume of the polyhedron P  using V-representation. If the
%  polyhedron is not in V-representation, it will be converted automatically.
%  Volume is inf if P  is unbounded, and zero if P  is not full-dimensional.
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
%          vol Volume of P.                             
%              Class: double                            
%                
%  
%  
%  SEE ALSO
%  --------
%     minVRep
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
 
 
no = numel(P);
if no>1
    vol = zeros(size(P));
    for i=1:no
        vol(i) = P(i).volume;
    end
    return;
end

% check emptyness and full-dimensionality
if P.isEmptySet || ~P.isFullDim
    vol = 0;
    return;
end

% computing volume requires vertices
P.minVRep();

if ~P.isBounded
	% unbounded polyhedra have infinite volume
	vol = Inf;
	
elseif P.Dim==1
	% issue #71: we need to handle 1D cases manually
	vol = max(P.V) - min(P.V);

elseif size(P.V, 1) == P.Dim+1
	% cheaper volume computation if the set is a fully-dimensional simplex
	% https://en.wikipedia.org/wiki/Simplex#Volume
	S = P.V';
	D = zeros(size(S, 1), size(S, 2)-1);
	for j = 2:size(S, 2)
		D(:, j-1) = S(:, j)-S(:, 1);
	end
	vol = 1/factorial(size(S, 2)-1)*abs(det(D));

else
	% general computation
	[~, vol] = convhulln(P.V);

end
