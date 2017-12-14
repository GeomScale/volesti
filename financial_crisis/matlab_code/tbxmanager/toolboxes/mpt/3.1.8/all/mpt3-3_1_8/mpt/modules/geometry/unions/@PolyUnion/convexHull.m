function H = convexHull(U)
%
%  CONVEXHULL: Computes the convex hull for union of polyhedra 
%  ============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      H = U.convexHull
%      H = convexHull(U)
%    
%  
%  DESCRIPTION
%  -----------
%     The convex hull of the union of polyhedra is defined as the minimal convex
%  set that contains all polyhedra. Note that computation of convex hull is an
%  expensive operation, therefore the result is stored internally under
%  Internal.convexHull which can be accessed.
%  
%  INPUT
%  -----
%     
%        
%          U Union of polyhedra in the same           
%            dimension.                               
%            Class: PolyUnion                         
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          H Convex hull of the polyhedra contained   
%            in the union                             
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     isConvex,  merge,  reduce
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

if numel(U)==0
	H = Polyhedron;
	return
	
elseif numel(U)>1
	% compute a single convex hull of all unions
	Hi = U.forEach(@(x) x.convexHull); % forEach() will use parfor
	H = PolyUnion(Hi).convexHull();
	return
end

% if there is 0 or 1 set contained, return
if U.Num==0
    H = Polyhedron;
    return
elseif U.Num==1
    H = Polyhedron(U.Set);
    return
end


if ~isfield(U.Internal,'convexHull') || isempty(U.Internal.convexHull)
    % compute the convex hull
	Vn=[];
	Rn=[];
	for i=1:U.Num
		if ~U.Set(i).hasVRep
			% if object is in H-rep, convert it to V
			U.Set(i).computeVRep();
		end
		Vn = [Vn; U.Set(i).V];
		Rn = [Rn; U.Set(i).R];
	end
	
	% construct one polyhedron
	H = Polyhedron('V',Vn,'R',Rn);
	
	% irredundant Vrep
	H.minVRep();
	
	% irredundant Hrep
	H.minHRep();
    
    % store internally
    U.Internal.convexHull = H;
else    
    H = U.Internal.convexHull;
end


end
