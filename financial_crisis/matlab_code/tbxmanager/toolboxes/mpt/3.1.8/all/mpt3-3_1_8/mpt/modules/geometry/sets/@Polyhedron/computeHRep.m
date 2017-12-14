function obj = computeHRep(obj)
%
%  COMPUTEHREP: Compute H-representation of a polyhedron. 
%  =======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      P = P.computeHRep
%    
%  
%  DESCRIPTION
%  -----------
%     Computes (possibly redundant) H-representation of the polyhedron: 
%                                         ---                   
%                    P = {x  |  Ax <= b } | | {x  |  A  x = b  }
%                                         | |         e      e  
%     This method implements facet enumeration using CDD solver and nlrs solver.
%  Please note that this is computationally demanding problem and the CDD solver
%  may become irresponsive for large input data.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in V-representation           
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          P Polyhedron in H-representation.          
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     computeVRep,  minVRep,  minHRep
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
 
 
global MPTOPTIONS

% deal with arrays
if numel(obj)>1
	if nargout==0
		obj.forEach(@computeHRep);
	else
		obj = obj.forEach(@computeHRep);
	end
	return
end

if obj.hasHRep
	% nothing to do
	return
elseif ~obj.hasVRep
	% empty set
	obj.hasHRep = true;
	return
elseif obj.isFullSpace()
	% R^n
	Rn = Polyhedron.fullSpace(obj.Dim);
	obj.H_int = Rn.H;
	obj.He_int = Rn.He;
	obj.hasHRep = true;
	return
end

% compute Hrep

obj.minVRep();
if isempty(obj.R_int) && obj.Dim>1 && obj.isFullDim() && ...
		size(obj.V_int, 1)>=obj.Dim+1 && size(obj.V_int, 2)<=3
	% try convhulln first; requires following conditions to be met:
	% * no rays
	% * dimension at least 2
    % * the set is full-dimensional
	% * the set has at least d+1 vertices

	x0 = obj.interiorPoint().x;
	V = obj.V_int;
    K = convhulln(V); % {'QJ', 'Qx', 'Qs'}
	d = size(V, 2);
	s.A = zeros(size(K, 1), d);
	s.B = ones(size(K, 1), 1);
	s.lin = [];
    for i = 1:size(K, 1)
        % each row of K contains indices of vertices that lie on the i-th
        % facet
        P = V(K(i, :), :);
        % compute the normal vector and the offset of the facet
        W = [P, -ones(d, 1)];
        [AB, ~] = qr(W'); % qr() is much faster than null()
        a = AB(1:d, end);
        b = AB(end, end);
        
        % determine the sign
        if a'*x0>b
            a = -a;
            b = -b;
        end
        s.A(i, :) = a';
        s.B(i) = b;
    end

else
    % Do facet enumeration with CDD
    s = cddmex('hull', struct('V', obj.V_int, 'R', obj.R_int));
end

Hall  = [s.A s.B];
H = Hall; H(s.lin, :) = [];
He = Hall(s.lin, :);
obj.H_int = H;
obj.He_int = He;
obj.hasHRep = true;

% unset obj.optMat since the H-representation might have changed
obj.optMat = [];

end
