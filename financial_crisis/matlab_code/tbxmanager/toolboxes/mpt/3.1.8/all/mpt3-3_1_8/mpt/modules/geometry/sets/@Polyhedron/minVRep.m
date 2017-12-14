function obj = minVRep(obj)
%
%  MINVREP: Compute an irredundant V-representation of a polyhedron. 
%  ==================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      P.minVRep()
%    
%  
%  DESCRIPTION
%  -----------
%     Computes an irredundant V-representation of the polyhedron.
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
%  SEE ALSO
%  --------
%     minHRep
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
 
 
if numel(obj)>1
	obj.forEach(@minVRep);
	return
end

if ~obj.hasVRep
	% automatically convert to Vrep if necessary
	obj.computeVRep();
end
if obj.irredundantVRep
	% nothing to do here
	return
end

if isempty(obj.V_int) && isempty(obj.R_int)
	% still not Vrep available, probably an empty set
	return
end

use_cddmex = true;
if isempty(obj.R_int) && obj.Dim>1 && obj.isFullDim() && ...
		size(obj.V_int, 1)>=obj.Dim+1
	% try convhulln first; requires following conditions to be met:
	% * no rays
	% * dimension at least 2
    % * the set is full-dimensional
	% * the set has at least d+1 vertices
    try
        K = convhulln(obj.V_int);
        s.V = obj.V_int(unique(K), :);
        s.R = obj.R_int;
        use_cddmex = false;
    end
end
if use_cddmex
    s = cddmex('reduce_v', struct('V', obj.V_int, 'R', obj.R_int));
end
	
obj.V_int = s.V;
obj.R_int = s.R;
obj.irredundantVRep = true;
% unset obj.optMat since the H-representation might have changed
obj.optMat = [];

% Lift the polyhedron to a cone
%R = [obj.R zeros(size(obj.R,1),1);obj.V ones(size(obj.V,1),1)];

% Detect if the affine hull is not Rn

%%% TODO : Do redundancy elimination via LP sequence

end
