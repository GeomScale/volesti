function alpha = shoot(obj, x)
%
%  SHOOT: Compute the maximal value of a multiplier in the desired direction. 
%  ===========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      alpha = S.shoot(x)
%      alpha = shoot(S, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the maximal value of the multiplier alpha  in the direction given by
%  the point x  by solving an optimization problem 
%                                         max   alpha            
%                              s.t.    alpha x in Set            
%     This problem is usually referred as shooting towards the point xbecause the
%  point alpha x  should lie on the edge of the set.
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
%          alpha The optimal value of the multiplier      
%                alpha such that alpha*x lies inside the  
%                Set. If the set is empty, not-a-number   
%                is returned. For unbounded direction, an 
%                inf-value is returned.                   
%                Class: double                            
%                  
%  
%  
%  SEE ALSO
%  --------
%     contains,  extreme,  project
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
no = numel(obj);
if no>1
    alpha = NaN(size(obj));
    for i=1:no
        alpha(i) = obj(i).shoot(x);        
    end
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

% do not put obj.alpha directly into constraints because YALMIP
% changed its status to quadratic scalar which caused later an error when 
% detecting appropriate solver
a = obj.alpha;
% solve the problem via YALMIP
F = obj.constraints + [obj.vars(:) == x(:)*a];
d = solvesdp(F, -a, obj.opts);

% if we don't know if it is feasible or unbounded, retry with artificial bounds
if ismember(d.problem,[12, 15])
    F = F + [ -MPTOPTIONS.infbound*ones(size(obj.vars)) <= obj.vars <= MPTOPTIONS.infbound*ones(size(obj.vars)) ];
    F = F + [ -MPTOPTIONS.infbound <= a <= MPTOPTIONS.infbound ];
    d = solvesdp(F, -a, obj.opts);
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
        alpha = double(a);
    case MPTOPTIONS.INFEASIBLE,
        alpha = NaN;
    case MPTOPTIONS.UNBOUNDED,
        alpha = Inf;
    otherwise
        error('Solver returned "%s" error when called from YALMIP.',sol.how);
end
end
