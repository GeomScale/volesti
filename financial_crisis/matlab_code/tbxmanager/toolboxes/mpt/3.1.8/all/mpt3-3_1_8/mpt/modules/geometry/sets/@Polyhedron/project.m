function sol = project(obj, y)
%
%  PROJECT: Project a point onto the given polyhedron. 
%  ====================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ret = P.project(x)
%      ret = project(P, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes the projection of the point y  onto this polyhedron by solving the
%  optimization problem: 
%                                         2             
%                             min{ ||y-x||   |  x in P }
%                                         2             
%     If P  is a vector of m  Polyhedra and x  is a matrix of size n X P.Dim, then
%  ret is a cell array with the structure fields that correspond the projection of
%  each point in the matrix.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          y Vector of size P.Dim that represent      
%            single point y, or matrix of size P.DimX 
%            nthat corresponds to multiple points     
%            merged column-wise.                      
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ret          Optimal solution                         
%                       Class: struct                            
%          ret.x        Projected point or [] if P  is empty     
%          ret.exitflag Informs about the termination status of  
%                       the optimization problem.                
%                       Class: double                            
%                       Allowed values:                          
%                                                                
%                         mptopt.OK                              
%                         mptopt.INFEASIBLE                      
%                         mptopt.UNBOUNDED                       
%                         mptopt.ERR                             
%                                                                
%          ret.dist     Distance from y  to the set P.           
%                         
%  
%  
%  SEE ALSO
%  --------
%     projection
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

% reject arrays
error(obj.rejectArray());

dim = obj.Dim;
if dim<1
    error('Cannot project with empty polyhedra.');
end
validate_realmatrix(y);
if ~isequal(size(y, 1), obj.Dim)
	error('Input argument must have %d rows.', obj.Dim);
end

%% Project points onto the polyhedra
n_points = size(y, 2);
sol(1, n_points) = struct('x',[],'exitflag',[],'dist',[]); 

for j = 1:n_points
    
    % (x-y)'(x-y) = x'x - 2*x'y + y'y
    qp   = obj.optMat;
 
    % semidefinite QP
    cost = obj.buildCost(-2*y(:,j), 2*eye(obj.Dim));
    qp.f = cost.f; qp.H = cost.H;    
	qp.quickqp = true;
    opt  = mpt_solve(qp);
    
    % if not feasible, retry with setting bounds on all variables - helps
    % to recover feasibility 
    if opt.exitflag ~= MPTOPTIONS.OK
       qp.lb(qp.lb<-MPTOPTIONS.infbound) =-MPTOPTIONS.infbound;
       qp.ub(qp.ub>MPTOPTIONS.infbound) = MPTOPTIONS.infbound;
       opt  = mpt_solve(qp);
    end
    
    sol(1,j).exitflag    = opt.exitflag;
    sol(1,j).dist    = Inf; % infeasible proble = infinite distance by convention
    if sol(1,j).exitflag == MPTOPTIONS.OK,
        sol(1,j).x       = opt.xopt(1:obj.Dim);
        sol(1,j).dist = norm(sol(1,j).x - y(:,j));
    end
end



end

