function tf = contains(P, x, fastbreak)
%
%  CONTAINS: Test if a polyhedron/point is contained inside polyhedron. 
%  =====================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      tf = P.contains(S)
%      tf = P.contains(S, fastbreak)
%      tf = contains(P, S, fastbreak)
%      P > S
%      P >= S
%      S < P
%      S <= P
%    
%  
%  DESCRIPTION
%  -----------
%     Check if the polyhedron or the point S is contained inside the polyhedron P.
%  The result it the logical statement if S subseteq P  and false otherwise. S can
%  be given as polyhedron, it can be a point given as a column (or set of points
%  concatenated column-wise). Notes: 
%    
%     - P.contains(S) with S a point assumes that S is a column vector. No
%     automatic transposition is performed! 
%     - P.contains(S) with S a matrix of points assumes that the points are stored
%     column-wise. No automatic transposition is performed! 
%     - If P is an array, P.contains(S) returns a (n_p X n_S) matrix of logicals,
%     where n_p is the number of elements of  P and n_x the number of points (i.e.,
%     the no. of columns) in S. 
%    Remember that testing the set membership is numerically sensitive and depends
%  on the settings of the absolute and relative tolerances.
%  
%  INPUT
%  -----
%     
%        
%          P         polyhedron or an array of polyhedra with 
%                    m elements                               
%                    Class: Polyhedron                        
%          S         Polyhedron or a polyhedron array or a    
%                    set of points. Multiple points must be   
%                    concatenated in a matrix column-wise,    
%                    where the number of rows corresponds to  
%                    the dimension of P and the number of     
%                    columns gives the number of points to    
%                    test. No automatic transposition is      
%                    performed!                               
%                    Class: Polyhedron or double              
%          fastbreak A logical flag. If true it indicates     
%                    that the set membership test should be   
%                    terminated as soon as at least one       
%                    element of P contains the point S.       
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      0                                      
%                      1                                      
%                                                             
%                    Default: 0                               
%                      
%  
%  
%  OUTPUT
%  ------
%     
%        
%          tf True if S subseteq P  and false          
%             otherwise.                               
%                                                      
%                - If S is a single point or a single  
%                polyhedron: status is a (m X 1)       
%                vector of true/false,  status(i)=true 
%                iff P(i) contains S.                  
%                - If S is a matrix of points (i.e., S 
%                = [s_1, ..., s_n]) status is a (m X   
%                n) matrix of true/false,  status(i,   
%                j)=true iff P(i) contains point S(:,  
%                j).                                   
%                - If S is a polyhedron, then          
%                P.contains(x) iff x subseteq P        
%                                                      
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
%     interiorPoint,  le,  lt,  ge,  gt
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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

narginchk(2, 3);
if ~( isnumeric(x) || isa(x, 'Polyhedron') )
	error('The input must be a real vector/matrix or a Polyhedron object.');
end
if nargin<3
	fastbreak = false;
end

m = numel(P);
if isnumeric(x)
	[d,n] = size(x);    
else
	n = numel(x);
	if n > 1
		% TODO: use setdiff for arrays
		error('Can only test containment of a single polyhedron.');
	end
end

if m==0
	% empty polyhedron array
	tf = [];
	
elseif isnumeric(x) && m==1 && P.hasHRep
	% special case:
	%   P = single polyhedron in H-rep
	%   x = single or multiple points
	%
	% this is a frequent case in Polyhedron/meshGrid, which needs to be as
	% fast as possible to have decent runtime of Polyhedron/fplot
    if d~=P.Dim
        error('The point(s) must be %dx1 vector(s).', P.Dim);
    end
	tf = false(m, n);
	for i = 1:n
		% iterate over points
		if any(P.H_int*[x(:, i); -1] > MPTOPTIONS.abs_tol)
			% not in the inequality Hrep
		elseif ~isempty(P.He_int) && ...
				any(abs(P.He_int*[x(:, i); -1]) > MPTOPTIONS.abs_tol)
			% not in the equality Hrep
		else
			tf(i) = true;
		end
	end
	
elseif isnumeric(x) && all([P.hasHRep])
	% special case:
	%   P = array of H-rep polyhedra
	%   x = single or multiple points
	%
	% this is so frequently used by PolyUnion/contains that it deserves a
	% fast implementation
	tf = false(m, n);
	for i = 1:n
		[~, inwhich] = P.isInside(x(:, i), struct('fastbreak', fastbreak));
		tf(inwhich, i) = true;
	end
	
elseif m>1
	% "P" is an array
	tf = false(m, n);
	for i = 1:m
		tf(i, :) = P(i).contains(x, fastbreak);
	end
	
elseif isnumeric(x) && n>1
	% "P" is a single polyhedron, "x" = multiple points
	tf = false(1, n);
	for i = 1:n
		tf(i) = P.contains(x(:, i), fastbreak);
	end
	
elseif isnumeric(x)
	% "P" is a single polyhedron, "x" is a single point
	tf = sub_contains_point(P, x);
	
else
	% "P" is a single polyhedron, "x" is a single polyhedron
	tf = sub_contains_polyhedron(P, x);
	
end

end

%-------------------------------------
function tf = sub_contains_polyhedron(P, S)

global MPTOPTIONS

if S.isEmptySet()
	% empty set is always contained in any other set
	tf = true;
	return
	
elseif P.isEmptySet()
	% empty set cannot contain a non-empty set (note that at this
	% point, due to the previous check of x.isEmptySet(), we know that
	% "x" is not empty)
	tf = false;
	return
	
elseif S.Dim ~= P.Dim
	error('Polyhedron S must be of the dimension %i.', P.Dim);
	
end

if isempty(S.R_int) && isempty(P.R_int) && ...
        ~(S.hasVRep && P.hasHRep)
	% check outer approximations first, but only if we don't have any rays,
	% otherwise computing the outer approximation is expensive. If S is in
	% V-rep and P in Hrep, we also bypass since P.contains(S) is cheap to
	% test directly.
	P.outerApprox();
	S.outerApprox();
	Plb = P.Internal.lb;
	Pub = P.Internal.ub;
	Slb = S.Internal.lb;
	Sub = S.Internal.ub;
	bbox_tol = MPTOPTIONS.rel_tol*1e2;
    if any(Slb + bbox_tol < Plb) || any(Sub - bbox_tol > Pub),
        % outer approximation of S is not contained in the outer
        % approximation of P, hence S cannot be contained in P
        tf = false;
        return
    end
end

tf = true;
if S.hasVRep && (P.hasHRep || (~isempty(S.V_int) && isempty(S.R_int)))
	% to check P.contains(S) with S in V-rep we reuquire either
	% 1) P to be in H-rep, or
	% 2) S to have vertices and no rays
	if P.hasHRep
		% Easiest case => test each ray and vertex
		P.minHRep();
		
		% check also equalities
		A = [P.A;P.Ae;-P.Ae]; b = [P.b;P.be;-P.be];
		
		% Let:
		%   P = { x | a'*x <= b }
		%   S = { x | x=r*y, y >= 0 }
		% Then P.contains(S) iff \forall x \in S we have x \in P, i.e.
		%   a'*r*y <= b, \forall y>=0
		% which holds iff
		%   a'*r <= min(0, b)
		
		I = A*S.R' - repmat(min(0, b), 1, size(S.R,1));
		if any(I(:) > MPTOPTIONS.rel_tol),
			% rays are not contained
			tf = false;
			return
		end
		
		I = A*S.V' - repmat(b,1,size(S.V,1));
		if any(I(:) > MPTOPTIONS.rel_tol),
			tf = false;
			return
		end
		
		% if all of the tests passed, check containment of a single point to
		% verify it is ok
		ip = S.interiorPoint;
		if ~P.contains(ip.x)
			tf = false;
			return
		end
		
	else
		% Test containment of each vertex of S in the points of P
		if any(~P.contains(S.V')) % points must be stored column-wise
			tf = false;
			return
		end
		
	end
else
	% S is an H-rep or in incompatible V-rep => Need an H-rep of P and S

    % only compute minimal H-representations after performing the
    % boundingbox-based heuristics
    P.minHRep();

	% check also equalities
	A = [P.A;P.Ae;-P.Ae]; b = [P.b;P.be;-P.be];
	
	% Easy case => support for each row of S must be less than that for P
	for i=1:size(A,1)
		s = S.extreme(A(i,:)');
		if s.exitflag~=MPTOPTIONS.OK
			% infeasible, or unbounded, or error  -> return false
			tf = false;
			return
		end
		% check also equalities if present
		if ~isempty(P.He_int)
			nhe = norm(P.Ae*s.x - P.be,Inf);
		else
			nhe=0;
		end
		
		% scaling the support is better if normalizing by b(i) -> see
		% test_polyhedron_contains_05_pass for P, 10*P, 100*P, etc.
		if norm(b(i),Inf)>MPTOPTIONS.rel_tol
			if ( s.supp/abs(b(i)) > b(i)/abs(b(i)) + MPTOPTIONS.rel_tol) || nhe>MPTOPTIONS.rel_tol
				tf = false;
				return
			end
		elseif ( s.supp > b(i) + MPTOPTIONS.rel_tol) || nhe>MPTOPTIONS.rel_tol
			tf = false;
			return
		end
	end
end

end

%-------------------------------------
function tf = sub_contains_point(P, x)

global MPTOPTIONS

nx = numel(x);
if P.isEmptySet()
	% empty set cannot contain any point
	tf = false;
	
elseif nx~=P.Dim
	error('The vector/matrix "x" must have %i rows.', P.Dim);

elseif P.hasHRep
	tf = P.isInside(x);
	
else
	% We only have V-rep, so do it the hard way
	sol = P.project(x);
	% normalized distance scales better for increasing the size of the
	% polyhedron -> try test_polyhedron_contains_04_pass for P, 10*P, 100*P, etc.
	tf = (sol.dist < MPTOPTIONS.rel_tol);

end

end
