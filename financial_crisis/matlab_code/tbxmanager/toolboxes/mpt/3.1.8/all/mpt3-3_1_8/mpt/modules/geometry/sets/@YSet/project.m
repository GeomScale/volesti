function sol = project(obj, x)
%
%  PROJECT: Compute the projection of the point onto this set. 
%  ============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      s = S.project(x)
%      s = project(S, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the projection of the point x onto this set. The projection problem
%  can be also explained as finding the closest point y from the set S to the given
%  point x.
%  
%  INPUT
%  -----
%     
%        
%          S A convex set described as YSet object.   
%            Class: YSet                              
%          x A point given as vector. Note that for   
%            YSet with symmetric matrix variable, the 
%            point x must be given as vector with     
%            symmetric terms.                         
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          s          The output structure with the            
%                     information about the projected point    
%                     and the exit status from the             
%                     optimization.                            
%                     Class: struct                            
%          s.exitflag Exit status from the optimization, i.e.  
%                     an integer value that informs if the     
%                     result was feasible (1), or otherwise    
%                     (different from 1).                      
%                     Class: double                            
%          s.how      A string that informs if the result was  
%                     feasible ('ok'), or if any problem       
%                     appeared through optimization.           
%                     Class: char                              
%          s.x        Projected point that lies inside the set 
%                     S.                                       
%                     Class: double                            
%          s.dist     The distance between the projected point 
%                     and the given point x.                   
%                     Class: double                            
%                       
%  
%  
%  SEE ALSO
%  --------
%     contains,  extreme,  shoot
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

narginchk(2, 2);

% check vector x
validate_realvector(x);

% reject arrays
error(obj.rejectArray());

% check dimension
if ~isequal(size(obj.vars), size(x))
	error('Input argument "x" must be a %dx%d matrix.', ...
		size(obj.vars, 1), size(obj.vars, 2));
end

% if x was created out of the matrix, there might be symmetric terms, 
% the assign command checks for the compatibility of the vector
assign(obj.vars, x);

% solve the problem via YALMIP
cost = norm(x - obj.vars,2)^2;
d = solvesdp(obj.constraints, cost, obj.opts);

% if we don't know if it is feasible or unbounded, retry with artificial bounds
if ismember(d.problem,[12, 15])
    F = obj.contraints + [ -MPTOPTIONS.infbound*ones(size(obj.vars)) <= obj.vars <= MPTOPTIONS.infbound*ones(size(obj.vars)) ];
    d = solvesdp(F, cost, obj.opts);
    % if solution is feasible -> unbounded
    if d.problem == 0
        d.problem = 2;
    else
        % infeasible
        d.problem = 1;
    end
end

% get MPT flags
sol = yalmip2mptflag(d);

switch sol.exitflag
    case MPTOPTIONS.OK,
        sol.x = double(obj.vars);
        sol.dist = sqrt(double(cost));
    case MPTOPTIONS.INFEASIBLE,
        sol.x = [];
        sol.dist = NaN;
    case MPTOPTIONS.UNBOUNDED,
        sol.x = [];
        sol.dist = Inf;
    otherwise
        error('Solver returned "%s" error when called from YALMIP.',sol.how);
end
end
