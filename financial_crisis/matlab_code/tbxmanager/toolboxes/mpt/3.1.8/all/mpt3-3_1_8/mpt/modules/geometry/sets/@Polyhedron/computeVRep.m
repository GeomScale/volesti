function obj = computeVRep(obj)
%
%  COMPUTEVREP: Compute V-representation of a polyhedron. 
%  =======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      P = P.computeVRep
%    
%  
%  DESCRIPTION
%  -----------
%     Computes (possibly redundant) V-representation of the polyhedron: 
%               P = V'lam | lam >= 0, sum(lam) = 1 + R'gam | gam >= 0 
%    This method implements vertex enumeration using CDD solver and nlrs solver.
%  Please note that this is computationally demanding problem and the CDD solver
%  may become irresponsive for large input data.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in H-representation           
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          P Polyhedron in V-representation.          
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     computeHRep,  minVRep,  minHRep
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
		obj.forEach(@computeVRep);
	else
		obj = obj.forEach(@computeVRep);
	end
	return
end

if obj.hasVRep
	% nothing to do
	return
elseif ~obj.hasHRep
	% empty set
	obj.hasVRep = true;
	return
elseif obj.isEmptySet()
	% empty set = empty vertices and rays
	obj.V_int = zeros(0, obj.Dim);
	obj.R_int = zeros(0, obj.Dim);
	obj.hasVRep = true;
	return
elseif obj.isFullSpace()
	% R^n = zero vertex and all basis vectors as rays
	obj.V_int = zeros(1, obj.Dim);
	obj.R_int = [eye(obj.Dim); -eye(obj.Dim)];
	obj.hasVRep = true;
	return
end

done = false;
backupTried = false;

% work with minimal H-representations to improve numerics
obj.minHRep();

% shift the polytope such that it contains the origin in its interior --
% helps CDD a lot
xc = zeros(obj.Dim, 1);
if obj.isBounded()
    xc = obj.chebyCenter.x;
end
Ai = obj.A;
bi = obj.b - Ai*xc;
Ae = obj.Ae;
be = obj.be - Ae*xc;
A = [Ae; Ai]; % equalities must come first
B = [be; bi];

% Do vertex enumeration with CDD
while ~done
	try
		% change almost-zero elements to zero
		A(abs(A)<MPTOPTIONS.zero_tol) = 0;
		B(abs(B)<MPTOPTIONS.zero_tol) = 0;
		
		% round H-representation to certain number of decimal places
		roundfactor = 10^15;
		A = round(A*roundfactor) / roundfactor;
		B = round(B*roundfactor) / roundfactor;
		
		lin = 1:size(obj.He_int, 1);
		% call CDD
		try
			s = cddmex('extreme', struct('A',A,'B',B,'lin',lin));
		catch
			% if CDD fails, retry with Matlab version
			[s.V,s.R,s.A,s.B] = mpt_nlrs('extreme',obj);
        end
        % shift vertices back
		s.V = s.V + repmat(xc', size(s.V, 1), 1);
        
		if size(s.V,1) == 0 % This is a cone... we need an explicit vertex
			obj.V_int = [obj.V_int; zeros(1,obj.Dim)];
		else
			obj.V_int = [obj.V_int; s.V];
		end
		obj.R_int = [obj.R_int; s.R];
		done = true;
	catch
		if backupTried
			obj.V_int = [obj.V_int; zeros(1,obj.Dim)];
			done=true;
			% this error appears usually when plotting polyhedra,
			% it is therefore disabled to show at least something
			% error('Could not compute vertices : Numerical problems in CDD.');
		end
		backupTried = true;
		
		% Use trick from mpt2.6 and reduce the calculation precision
		% CDD sometimes fails to compute extreme points correctly
		% this happens often when slopes of two hyperplanes are too close
		% that's why we use a fixed-point arithmetics
		
		thull = Polyhedron('H', fix(obj.H_int * 1e5) / 1e5, ...
			'He', fix(obj.He_int*1e5) / 1e5);
		thull.minHRep();
		% Do vertex enumeration with CDD
		A = [thull.He_int(:,1:end-1); thull.H_int(:,1:end-1)];
		B = [thull.He_int(:,end); thull.H_int(:,end)];
	end
end

obj.hasVRep = true;

% unset obj.optMat since the H-representation might have changed
obj.optMat = [];

end
