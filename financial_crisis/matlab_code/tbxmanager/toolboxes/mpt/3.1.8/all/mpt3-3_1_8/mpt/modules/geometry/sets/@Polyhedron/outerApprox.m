function approx = outerApprox(obj)
%
%  OUTERAPPROX: Computes outer bounding hypercube of a polyhedron. 
%  ================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      B = P.outerApprox
%      P.outerApprox
%    
%  
%  DESCRIPTION
%  -----------
%     B = P.outerApprox computes the smallest axis-aligned hypercube B that
%  contains the polyhedron P. The lower and upper bounds of the hypercube are also
%  stored in P.Internal.lb and P.Internal.ub, respectively.
%    Use P.outerApprox if you only want the bounds to be written to P.Internal, but
%  do not need the explicit bounding hypercube to be constructed.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron.                              
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          B Bounding hypercube.                      
%            Class: Polyhedron                        
%              
%  
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

if length(obj)>1
	for i = 1:length(obj)
		if nargout==1
			approx(i) = obj(i).outerApprox;
		else
			obj(i).outerApprox;
		end
	end
	return
end

d = obj.Dim;

if isfield(obj.Internal, 'lb') && isfield(obj.Internal, 'ub')
	% reuse stored information
	lb = obj.Internal.lb;
	ub = obj.Internal.ub;

elseif obj.hasVRep
	% Vrep is easy without rays, just take min/max of vertices
	if isempty(obj.R_int)
		lb = min(obj.V_int, [], 1)';
		ub = max(obj.V_int, [], 1)';
	else
		% resort to ConvexSet/outerApprox for unbounded polyhedra
		approx = outerApprox@ConvexSet(obj);
		lb = approx.Internal.lb;
		ub = approx.Internal.ub;
	end
	% update properties of the input object
	obj.Internal.lb = lb;
	obj.Internal.ub = ub;
	
elseif obj.isEmptySet()
    lb = Inf(obj.Dim, 1);
    ub = -Inf(obj.Dim, 1);
    obj.Internal.lb = lb;
    obj.Internal.ub = ub;
    
else
	% for Hrep we have to solve 2 LPs per each dimension
	H = obj.H_int;
	He = obj.He_int;
	f = zeros(1, d);
	LP.f = f;
	LP.A = H(:, 1:end-1);
	LP.b = H(:, end);
	LP.Ae = He(:, 1:end-1);
	LP.be = He(:, end);
	LP.lb = []; 
	LP.ub = [];
	LP.quicklp = true;
	
	lb = -Inf(d, 1);
	ub = Inf(d, 1);
	for i = 1:d
		% minimize
		LP.f = f;
		LP.f(i) = 1;
		sol = mpt_solve(LP);
		if sol.exitflag == MPTOPTIONS.OK
			lb(i) = sol.obj;
		end
		
		% maximize
		LP.f(i) = -1;
		sol = mpt_solve(LP);
		if sol.exitflag == MPTOPTIONS.OK
			ub(i) = -sol.obj;
		end
	end
	
	% update properties of the input object
	obj.Internal.lb = lb;
	obj.Internal.ub = ub;
end

% construct output arguments
if nargout==1
	% return the bounding hyperrectangle as a Polyhedron object
	approx = Polyhedron([eye(d); -eye(d)], [ub; -lb]);
	approx.irredundantHRep = true;
	approx.Internal.lb = lb;
	approx.Internal.ub = ub;
end
