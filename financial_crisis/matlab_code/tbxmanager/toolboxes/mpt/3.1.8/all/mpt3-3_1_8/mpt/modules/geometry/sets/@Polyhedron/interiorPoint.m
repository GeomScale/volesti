function sol = interiorPoint(P, facetI)
%
%  INTERIORPOINT: Compute a point in the relative interior of the Polyhedron. 
%  ===========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      sol = P.interiorPoint
%      sol = P.interiorPoint(facetIndex)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute a point in the relative interior of the polyhedron. If facetIndex in
%  {1,...,P.size(P.A,1)}  is specified, then a point in the relative interior of P 
%  {x |  P.A(facetIndex,:)x <= P.b(facetIndex)}  is returned.
%  
%  INPUT
%  -----
%     
%        
%          P          Polyhedron in any format                 
%                     Class: Polyhedron                        
%          facetIndex Index of an inequality of P  (row of     
%                     P.H).                                    
%                     Class: integer                           
%                       
%  
%  
%  OUTPUT
%  ------
%     
%        
%          sol                                                   
%                       Class: struct                            
%          sol.x        The interior point                       
%                       Class: double vector                     
%          sol.isStrict The output is true if x  is in the       
%                       strict relative interior, false          
%                       otherwise.                               
%                       Class: logical                           
%          sol.r        Radius of the largest ball centered at x 
%                        that is still within P                  
%                            y - x in P, forall ||y|| <= r       
%                        Note : r  is empty if P  is empty or    
%                       only has a V-rep.                        
%                       Class: double                            
%                         
%  
%  
%  SEE ALSO
%  --------
%     chebyCenter,  facetInteriorPoints
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

if nargin < 2
    facetI = [];
end

%% allocate output
% must be here to always have the same ordering of fields
sol = struct('x', [], 'isStrict', [], 'r', []);

%% deal with arrays
if numel(P)>1
	% return an array of structures
	sol = P.forEach(@(elem) elem.interiorPoint(facetI));
    return
end

%% try to reuse a stored information
if isempty(facetI) && ~isempty(P.Internal.InternalPoint)
	sol = P.Internal.InternalPoint;
	return
end	

%% checks for polyhedra

% if facets provided, P must be irredundant
if ~isempty(facetI)
	if ~P.irredundantHRep
		error('Polyhedron must be in minimal representation when you want compute any interior point of its facets. Use "minHRep()" to compute the facets.');
	end
end

%% compute interior points

% try specific problems first for faster computation
if isempty(facetI)
    if P.hasVRep
        % Average vertices and rays
		if isempty(P.V) && ~isempty(P.R)
			% only rays
			sol.x = mean(normalize(P.R), 1)';
		else
			% vertices and maybe rays
			sol.x = mean(P.V,1)';
			if size(P.R,1) > 0
				sol.x = sol.x + norm(sol.x)*mean(normalize(P.R), 1)';
			end
		end
        
        V = P.V;
		if ~isempty(V)
			V = V(2:end,:) - repmat(V(1,:),size(V,1)-1,1);
		end
        R = P.R;
        
        sol.isStrict = false;
        if rank([V;R], MPTOPTIONS.abs_tol) == P.Dim
            sol.isStrict = true;
		end
		if isempty(V)
			% the set is unbounded
			sol.r = Inf;
		end
		P.Internal.InternalPoint = sol;
        return
    end
    
    % It's an affine set
    if size(P.He_int,1) > 0 && isempty(P.H_int)
        % only feasible sets
        if size(P.He,1)<=P.Dim
            sol.x = P.Ae\P.be;
            sol.isStrict = false;
            sol.r = inf;
			P.Internal.InternalPoint = sol;
			return
        end
    end
else
	% check if facet is vector of indices
	validate_indexset(facetI);
end

% must test for H-rep in case P is empty
if P.hasHRep
    % by default call chebyCenter for all problems
    chb = P.chebyCenter(facetI);
    if chb.exitflag == MPTOPTIONS.OK
        sol.x = chb.x;
        sol.r = chb.r;
        if 2*chb.r > MPTOPTIONS.region_tol
            sol.isStrict = true;
        else
            sol.isStrict = false;
        end
        if ~isempty(P.He_int)
            sol.isStrict = false;
		end
		if isempty(facetI)
			P.Internal.InternalPoint = sol;
		end
    end
end


end
