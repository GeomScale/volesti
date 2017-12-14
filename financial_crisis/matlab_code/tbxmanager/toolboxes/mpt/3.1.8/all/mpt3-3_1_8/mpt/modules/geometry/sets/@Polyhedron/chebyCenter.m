function sol = chebyCenter(P, facetI, bound)
%
%  CHEBYCENTER: Compute the Chebyshev centre of the Polyhedron. 
%  =============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      s = P.chebyCenter
%      s = P.chebyCenter(facet)
%      s = P.chebyCenter(facet,bound)
%      s = chebyCenter(P,facet,bound)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes the Chebyshev centre of a polyhedron via solving the following
%  optimization problem 
%                            x,r  max   r                                  
%                              T                                           
%                    s.t.    a  x + ||a || r <= b ,    i=1,...,m           
%                             i        i  2      i                         
%                                                                          
%                             A x = b                                      
%                              e     e                                     
%     where x  is the centre of the ball inscribed inside the H-representation of
%  the polyhedron Pand r  is the radius of the ball. The polyhedron is given as P =
%  { x  |  Ax <= b, A_ex=b_e }  with m  inequalities and m_e  equalities. It is
%  possible to compute the Chebyshev centre restricted to one or more facets of
%  polyhedron by providing vector of indices facet. Note that the polyhedron must
%  be in its minimal representation in order to have facets. If the facet indices
%  are provided, the inequalities in the above optimization problem belonging to
%  these indices change to equality constraints. It is also possible to upper bound
%  the radius r by specifying bound input. In this case the constraint  r<=bound 
%  is added.
%  
%  INPUT
%  -----
%     
%        
%          P     Polyhedron in H-representation.          
%                Class: Polyhedron                        
%          facet Vector of indices specifying a subset of 
%                facets for computing the Chebyshev       
%                centre.                                  
%                Class: double                            
%                Default: []                              
%          bound Upper bound on the radius of the         
%                inscribed ball.                          
%                Class: double                            
%                Default: Inf                             
%                  
%  
%  
%  OUTPUT
%  ------
%     
%        
%          s          Structure with outputs from the          
%                     Chebyshev centre optimization problem.   
%                     Class: struct                            
%          s.exitflag Informs if the optimization problem was  
%                     terminated properly.                     
%                     Class: double                            
%          s.x        Centre of the ball inscribed in          
%                     polyhedron P.                            
%                     Class: double                            
%          s.r        Radius of the ball inscribed in          
%                     polyhedron P.                            
%                     Class: double                            
%                       
%  
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

% Artificial bound on the radius. Logic:
% 1) Solve the LP with this big bound on the radius.
% 2) If the LP is infeasible => exit, infeasible.
% 3) If the radius is smaller than big_bound => exit, bounded polytope.
% 4) If the radius is larger, we either have a larger polytope, or an
%    unbounded polyhedron. Drop the bound altogether and solve the LP
%    again. Then:
% 4.1) If the LP without bounds on the radius is feasible => polytope
% 4.2) Otherwise => unbounded polyhedron
%
% Note: this logic is skipped if user specifies his own non-inf bound
BIG_BOUND = 1.234e6;

if nargin < 2,
	facetI = [];
	bound = Inf;
elseif ~isempty(facetI)
	% check if facet is vector of indices
	validate_indexset(facetI);
end
if nargin < 3,
    bound = Inf;
else
    if ~isscalar(bound)
        error('Upper bound the the radius must be scalar.');
    end
    if ~isnumeric(bound) || bound<=0 || ~isreal(bound) || isnan(bound)
       error('Upper bound on the radius must be a numeric real scalar, positive valued.');
    end
end

%% deal with arrays
if numel(P)>1
	% return an array of structures
	sol = P.forEach(@(elem) elem.chebyCenter(facetI, bound));
    return
end

%% allocate output
% must be here to always have the same ordering of fields
sol = struct('exitflag', [], 'x', [], 'r', -inf);

%% checks for polyhedra
if ~P.hasHRep
	if ~P.hasVRep
		% empty polyhedron
		sol.exitflag = MPTOPTIONS.INFEASIBLE;
		sol.x = [];
		sol.r = -inf;
		return;
	end
	% get the Hrep
	P.minHRep();
end

H = P.H_int; He = P.He_int;

% if facets provided, P must be irredundant
if ~isempty(facetI)
    if ~P.irredundantHRep
        error('Polyhedron must be in minimal representation when you want compute Chebyshev centre any of its facets. Use "minHRep()" to compute the facets.');
    end
end

%% if Chebyshev data are stored internally, use this information
if isfield(P.Internal,'ChebyData') && isempty(facetI) && isinf(bound)
    
    % only available if there is no facet and bound provided
    if ~isempty(P.Internal.ChebyData)
        sol = P.Internal.ChebyData;
        return
    end
end

% Chebyshev centre solves a following LP:
%
%    max r
% s.t.:   a_i'*xc + ||a_i||_2*r <= b_i    i=1,...,m
%         Ae*xc = be
%
% where a_i is the i-th row of the inequality matrix A in A*xc<=b and
% ||a_i||_2 is the 2-norm of that row


% inequality constraints on [xc, r]
% [ a_i, ||a_i||_2 ], i=1,...,m
A = H(:, 1:end-1);
S.A  = [A sqrt(sum(A.*A,2))];
S.b = H(:, end);

% equality constraints
% [ a_i   0 ], i=1, ..., me
Ae = He(:, 1:end-1);
S.Ae = [Ae zeros(size(Ae,1),1)];
S.be = He(:, end);

% if we want to compute center of the facet, add this facet to equality
% constraints
if ~isempty(facetI)
    if any(facetI > size(S.A,1))
        error('Facet index must be less than number of inequalities (%i).', size(S.A,1));
    end
    S.Ae = [S.Ae;S.A(facetI,1:end-1) zeros(length(facetI),1)];
    S.be = [S.be;S.b(facetI)];
    S.A(facetI,:) = [];
    S.b(facetI,:) = [];
end

% lower bounds on [xc, r]
S.A = [S.A; zeros(1,P.Dim), -1];
S.b = [S.b; 0];
% upper bounds
inf_bound = isinf(bound);
if inf_bound
	% always introduce at least an artificial upper bound on the radius
	bound = BIG_BOUND;
end
% make sure the bound on the radius is always the last constraint!
S.A = [S.A; zeros(1,P.Dim), 1];
S.b = [S.b; bound];

% the last value is -1 because it is maximization
S.f = [zeros(size(H, 2)-1,1);-1];

% solve problem
S.lb = []; S.ub = []; S.quicklp = true;
ret = mpt_solve(S);

% set output variables
sol.exitflag = ret.exitflag;

if ret.exitflag == MPTOPTIONS.OK
	if -ret.obj>MPTOPTIONS.zero_tol
		sol.x = ret.xopt(1:end-1);
		sol.r = ret.xopt(end);
	else
		sol.x = ret.xopt(1:end-1);
		sol.r = 0;
	end
	if inf_bound && sol.r >= BIG_BOUND-1
		% The radius is too large, potentially we have an unbounded
		% polyhedron. Re-solve the LP without bounds to be sure.
		%
		% Introduce a numerical tolerance when comparing the radius to
		% BIG_BOUND.
		%
		% Remove the bound. If we get a bounded solution, then we have a
		% polytope. Otherwise we have a polyhedron.
		%
		% NOTE: as noted above, make sure the bounds are always the last
		% constraint!
		S.A = S.A(1:end-1, :);
		S.b = S.b(1:end-1);
		
		% Note that the LP can also be infeasible. But that indicates
		% unboundedness, since it was feasible with the bound included.
		ret = mpt_solve(S);
		if ret.exitflag == MPTOPTIONS.OK
			% bounded solution even without bounds on the radius => polytope
			sol.x = ret.xopt(1:end-1);
			sol.r = ret.xopt(end);
		else
			% unbounded => polyhedron
			sol.r = Inf;
		end
    end
    % check if the point is contained inside the polyhedron
    if ~P.contains(sol.x)
        sol.exitflag = MPTOPTIONS.INFEASIBLE;
    end
    
else
	% infeasible
	sol.x = [];
	sol.r = -Inf;
end

% store internally if empty facet and the bound was not provided
if isempty(facetI) && inf_bound
    P.Internal.ChebyData = sol;
end



end
