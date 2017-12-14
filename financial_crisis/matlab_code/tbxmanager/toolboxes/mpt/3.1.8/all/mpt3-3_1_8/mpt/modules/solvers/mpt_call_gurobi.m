function R = mpt_call_gurobi(S)
%
%  MPT_CALL_GUROBI: A gateway function to GUROBI solver (without errorchecks) 
%  ===========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_gurobi(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements call to GUROBI solver based on formulation from Opt
%  class. QP, LP, MILP and MIQP problems are supported. It is assumed that
%  QP/LP/MIQP/MILP entering this function (for LP/MILP H=0) is of the form 
%                                   1  T    T                      
%                             min   - x Hx+f x                 (1) 
%                                   2                              
%                            s.t.   lb <= x <= ub              (2) 
%                                   Ax <= b                    (3) 
%                                                                  
%                                   A x = b                    (4) 
%                                    e     e                       
%                                   x in {C, I, B, N, S }      (5) 
%     where the set {C, I, B, N, S}  represents 
%    
%     - C - continuous variables, x in (-oo,oo)   
%     - I - integer variables x in (..., -1, 0, 1, ...)   
%     - B - binary variables x in {0,1}  
%     - N - semi-integer variables (possibly bounded above)  x in [0, 1, x <= oo ) 
%     
%     - S - semi-continuous variables (possibly bounded above)  x in [0, x <= oo)  
%    which is given by strings in vartype field. GUROBI accepts this format
%  directly, the only prerequisite is to transform input data to sparse format.
%  
%  INPUT
%  -----
%     
%        
%          S              Structure of the Opt class.              
%                         Class: struct                            
%          S.H            Quadratic part of the objective          
%                         function.                                
%                         Class: double                            
%                         Default: []                              
%          S.f            Linear part of the objective function.   
%                         Class: double                            
%          S.A            Linear part of the inequality            
%                         constraints Ax <= b.                     
%                         Class: double                            
%          S.b            Right hand side of the inequality        
%                         constraints Ax <= b.                     
%                         Class: double                            
%          S.Ae           Linear part of the equality constraints  
%                         A_ex=b_e.                                
%                         Class: double                            
%                         Default: []                              
%          S.be           Right hand side of the equality          
%                         constraints A_ex=b_e.                    
%                         Class: double                            
%                         Default: []                              
%          S.lb           Lower bound for the variables x >= lb.   
%                         Class: double                            
%                         Default: []                              
%          S.ub           Upper bound for the variables x <= ub.   
%                         Class: double                            
%                         Default: []                              
%          S.n            Problem dimension (number of variables). 
%                                                                  
%                         Class: double                            
%          S.m            Number of inequalities in Ax <= b.       
%                         Class: double                            
%          S.me           Number of equalities in A_ex=b_e         
%                         Class: double                            
%          S.problem_type A string specifying the problem to be    
%                         solved.                                  
%                         Class: char                              
%          S.vartype      A string specifying the type of          
%                         variable. Supported characters are C     
%                         (continuous), I (integer), B (binary), N 
%                         (semi-integer), S (semi-continuous).     
%                         Example: First variable from three is    
%                         binary, the rest is continuous:          
%                         S.vartype='BCC';                         
%                         Class: char                              
%          S.test         Call (false) or not to call (true) MPT   
%                         global settings.                         
%                         Class: logical                           
%                         Default: false                           
%                           
%  
%  
%  OUTPUT
%  ------
%     
%        
%          R          result structure                         
%                     Class: struct                            
%          R.xopt     optimal solution                         
%                     Class: double                            
%          R.obj      Optimal objective value.                 
%                     Class: double                            
%          R.lambda   Lagrangian multipliers.                  
%                     Class: double                            
%          R.exitflag An integer value that informs if the     
%                     result was feasible (1), or otherwise    
%                     (different from 1).                      
%                     Class: double                            
%          R.how      A string that informs if the result was  
%                     feasible ('ok'), or if any problem       
%                     appeared through optimization.           
%                     Class: char                              
%                       
%  
%  
%  SEE ALSO
%  --------
%     mpt_solve
%  

%  AUTHOR(s)
%  ---------
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
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

assert(~isequal(S.problem_type, 'LCP'), 'mpt_call_gurobi: GUROBI solver does not solve LCP problems!');

% A*x <= rhs
model.A = sparse([S.Ae; S.A]);
model.rhs = full([S.be; S.b]);
model.sense = char(['='*ones(S.me, 1); '<'*ones(S.m, 1)]); 

% lb <= x <= ub
if isempty(S.lb)
    model.lb = -Inf(S.n, 1);
else
    model.lb = S.lb;
end
if isempty(S.ub)
    model.ub = Inf(S.n, 1);
else
    model.ub = S.ub;
end

% types of variables (Binary, Integer, ...)
if isfield(S, 'vartype') && ~isempty(S.vartype)
    model.vtype = S.vartype;
end

% minimize x'*H*x + obj'*x
model.obj = S.f;
if isequal(S.problem_type(end-1:end), 'QP')
    % quadratic term for QP/MIQP problems
    model.Q = sparse(S.H*0.5);
end

% set options
if S.test
    opts.OutputFlag = 0; % no verbosity
else
    opts = MPTOPTIONS.modules.solvers.gurobi;
end
% prevent gurobi from returning the INF_OR_UNBD status:
% http://www.gurobi.com/documentation/5.6/reference-manual/infunbdinfo
opts.InfUnbdInfo = 1;

% solve
result = gurobi(model,opts);

switch result.status
    case 'OPTIMAL'
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
        
        % objective value
        R.obj = result.objval;
        
        % primal optimizer
        R.xopt = result.x;
        
        % dual optimizer
        if ~isfield(result, 'pi')
            % no dual variables for mixed-integer problems
            lambda = NaN(size(model.A, 1), 1);
        else
            lambda = -result.pi;
        end
        
    case 'UNBOUNDED'
        R.how = 'unbounded';
        if S.test
            R.exitflag = 3;
        else
            R.exitflag = MPTOPTIONS.UNBOUNDED;
        end
        R.obj = -Inf; % unbounded minimization
        R.xopt = NaN(S.n, 1);
        lambda = NaN(size(model.A, 1), 1);
        
    case {'INFEASIBLE', 'INF_OR_UNBD'}
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end
        R.obj = Inf; % infeasible minimization
        R.xopt = NaN(S.n, 1);
        lambda = NaN(size(model.A, 1), 1);
        
    otherwise
        R.how = result.status;
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
        R.obj = Inf;
        R.xopt = NaN(S.n, 1);
        lambda = NaN(size(model.A, 1), 1);
        
end

% dual variables
R.lambda.ineqlin = lambda(S.me+1:S.me+S.m);
R.lambda.eqlin = lambda(1:S.me);
R.lambda.lower = NaN(S.n, 1);
R.lambda.upper = NaN(S.n, 1);

end
