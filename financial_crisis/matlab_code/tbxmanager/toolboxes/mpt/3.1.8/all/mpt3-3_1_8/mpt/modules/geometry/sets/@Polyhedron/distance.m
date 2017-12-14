function ret = distance(P, S)
%
%  DISTANCE: Compute the distance between the given point/polyhedron and this
%  ==========================================================================
%  polyhedron. 
%  ============
%  
%  
%  SYNTAX
%  ------
%     
%      dist = P.distance(x)
%      ret = P.distance(S)
%      dist = distance(P, x)
%      ret = distance(P, S)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the distance between the polyhedron P and the point x or the
%  polyhedron S. 
%    
%     1. 		By providing real vector x, the distance between x and P is computed by
%     		solving the optimization problem 		 			
%                                           2           
%                                min{||x-y||  | y in P }
%                                           2           
%    		 and the distance is returned as real number. 
%     2. 		If polyhedron S is specified as the argument, the distance between S and
%     P is 		computed by solving the following optimization problem 		 			
%                                       2                  
%                            min{||x-y||  | y in P, x in S}
%                                       2                  
%    		 where the results of the optimization are returned in a struct format. 
%  
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          x Vector of size P.Dim                     
%            Class: double vector                     
%          S Polyhedron with the same dimension as P. 
%                                                     
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          dist         Distance between the point x  and the    
%                       Polyhedron P                             
%                       Class: double                            
%          ret          Optimal solution or [] if P  is empty.   
%                       Class: struct                            
%          ret.exitflag Integer value informing about the        
%                       termination status of the optimization.  
%                       Class: double                            
%          ret.dist     Distance from S  to the set P            
%                       Class: double                            
%          ret.x        Point x  in S  closest to P.             
%                       Class: double                            
%          ret.y        Point y  in P  closest to S.             
%                       Class: double                            
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

assert(isa(P, 'Polyhedron'), 'The first input must be a Polyhedron object.');

% check if S is an array or polyhedron
if isa(S,'Polyhedron')
    if numel(S)>1
       error('Only single polyhedron S allowed.'); 
    end
else
    validate_realvector(S);
end

if numel(P)>1
	% return an array of structures
	ret = P.forEach(@(elem) elem.distance(S));
    return
end

% pre-alocate output
ret = struct('exitflag', [], 'dist', Inf, 'x', [], 'y', []);

%% S is supposedly a point
if ~isa(S,'Polyhedron') 
  % Call the superclass, which handles this case
  ret = P.distance@ConvexSet(S);
  return
end

%% S is Polyhedron
if S.Dim ~= P.Dim
    error('Both polyhedra have to be of the same dimension.');
end

if P.isEmptySet() || S.isEmptySet()
    % distance from an empty set is infinite by convention (issue #111)
    ret.exitflag = MPTOPTIONS.INFEASIBLE;
    ret.dist = Inf;
    return
end

% Get representations of both polyhedra
matP = P.optMat;
matS = S.optMat;

% Build optimization matrices
qp.A  = blkdiag(matP.A, matS.A);
qp.b  = [matP.b;matS.b];
qp.Ae = blkdiag(matP.Ae, matS.Ae);
qp.be = [matP.be; matS.be];
qp.lb = [matP.lb; matS.lb];
qp.ub = [matP.ub; matS.ub];

% Build the cost function (x-y)'(x-y)
nLamP = size(matP.A,2) - P.Dim;
nLamS = size(matS.A,2) - P.Dim;
HP  = diag([ones(P.Dim,1);zeros(nLamP,1)]);
HPS = blkdiag(-eye(P.Dim), zeros(nLamP, nLamS));
HS  = diag([ones(P.Dim,1);zeros(nLamS,1)]);
H   = [HP HPS; HPS' HS];
qp.H = H;
qp.f = [];

sol = mpt_solve(qp);
ret.exitflag = sol.exitflag;

if sol.exitflag==MPTOPTIONS.OK
    y = sol.xopt(1:P.Dim);
    x = sol.xopt(size(matP.A,2)+(1:P.Dim));
    ret.dist = norm(x-y);
    ret.x = x;
    ret.y = y;
end

end
