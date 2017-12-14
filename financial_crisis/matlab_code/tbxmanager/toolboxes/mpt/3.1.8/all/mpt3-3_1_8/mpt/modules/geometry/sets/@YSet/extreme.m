function sol = extreme(obj, x)
%
%  EXTREME: Compute an extreme point of this set in the given direction. 
%  ======================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      s = S.extreme(x)
%      s = extreme(S, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute an extreme point of this set in the direction given by the point x.
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
%                     information about the extreme point and  
%                     the exit status from the optimization.   
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
%          s.x        Computed extreme point that lies on the  
%                     boundary of the set S.                   
%                     Class: double                            
%          s.supp     The support of this set in the direction 
%                     x which represents the optimal value of  
%                     the objective function in the            
%                     optimization problem ymax   x^Ty,        
%                     s.t.    y in Set .                       
%                     Class: double                            
%                       
%  
%  
%  SEE ALSO
%  --------
%     contains,  project,  shoot
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

% deal with arrays
if numel(obj)>1
	% return an array of structures
	sol = obj.forEach(@(elem) elem.extreme(x));
    return
end

% check dimension
if numel(x)~=obj.Dim
    error('The argument must have %i number of elements.', obj.Dim);
end

if any(size(x)~=size(obj.vars))
    x = transpose(x);
end

% if x was created out of the matrix, there might be symmetric terms, 
% the assign command checks for the compatibility of the vector
assign(obj.vars, x);

% make column vector out of any matrix
x = x(:);

model = obj.extr.model;
model.c = 0*model.c;
model.c(obj.extr.local) = -x;
s  = feval(model.solver.call,model);

% if we don't know if it is feasible or unbounded, retry with artificial bounds
if ismember(s.problem,[12, 15])
    model.lb = -MPTOPTIONS.infbound*ones(size(model.lb));
    model.ub = MPTOPTIONS.infbound*ones(size(model.ub));
    s = feval(model.solver.call,model);
    
    % if solution is feasible -> unbounded
    if s.problem == 0
        s.problem = 2;
    else
        % otherwise infeasible
        s.problem = 1;
    end

end

% get MPT flags
sol = yalmip2mptflag(s);

switch sol.exitflag
    case MPTOPTIONS.OK
        v = s.Primal;
        v = v(obj.extr.local(:));
        sol.x       = v;
        sol.supp    = x'*v;
        if sol.supp >= MPTOPTIONS.infbound
            sol.exitflag = MPTOPTIONS.UNBOUNDED;
            sol.supp = Inf;
        elseif sol.supp <= -MPTOPTIONS.infbound
            sol.exitflag = MPTOPTIONS.UNBOUNDED;
            sol.supp = -Inf;            
        end
    case MPTOPTIONS.INFEASIBLE
        sol.x = [];
        sol.supp = NaN;
    case MPTOPTIONS.UNBOUNDED;
        sol.x = s.Primal(obj.extr.local(:));
        sol.supp = Inf;
    otherwise
        error('Solver returned "%s" error when called from YALMIP.',sol.how);
end

end
