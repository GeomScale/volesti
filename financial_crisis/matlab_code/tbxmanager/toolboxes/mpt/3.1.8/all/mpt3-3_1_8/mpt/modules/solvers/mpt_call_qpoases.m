function R = mpt_call_qpoases(S)
%
%  MPT_CALL_QPOASES: A gateway function to QPoases solver (without errorchecks) 
%  =============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_qpoases(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements call to QPoases solver based on formulation from Opt
%  class. Only QP and LP problems are supported. It is assumed that QP/LP entering
%  this function (for LP H=0) is of the form 
%                                      1  T    T                
%                                min   - x Hx+f x           (1) 
%                                      2                        
%                               s.t.   lb <= x <= ub        (2) 
%                                      Ax <= b              (3) 
%                                                               
%                                      A x = b              (4) 
%                                       e     e                 
%     which must be transformed to 
%                                     1  T    T                  
%                               min   - x Hx+f x             (5) 
%                                     2                          
%                              s.t.   lb <= x <= ub          (6) 
%                                                                
%                                     lA <= A x <= uA        (7) 
%                                            m                   
%                                                            (8) 
%     which accepts QPspline. Inequality (??)  and equality (??)  constraints are
%  merged to 
%                                ( -MPTOPTIONS.infbound  )
%                                (                       )
%                           lA = (          b            )
%                                (           e           )
%     
%                                         ( A   )
%                                         (     )
%                                    A  = ( A   )
%                                     m   (  e  )
%     
%                                         ( b   )
%                                         (     )
%                                    uA = ( b   )
%                                         (  e  )
%     since QPoases accepts equalities as double-sided inequalities.
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
%          S.me           Number of equalities in A_ex=b_e.        
%                         Class: double                            
%          S.problem_type A string specifying the problem to be    
%                         solved.                                  
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
%          R.xopt     Optimal solution.                        
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

if ~any(strcmpi(S.problem_type,{'QP','LP'}))
    error('mpt_call_qpoases: qpOASES solver does not solve %s problems!',S.problem_type);
end

% for LP we set H=0
if isempty(S.H)
    S.H = zeros(S.n);
end

% convert sparse matrices to full matrices
if issparse(S.H),
    S.H = full(S.H);
end
if issparse(S.f)
    S.f = full(S.f);
end
if issparse(S.A)
    S.A = full(S.A);
end
if issparse(S.b)
    S.b = full(S.b);
end
if issparse(S.Ae)
    S.Ae = full(S.Ae);
end
if issparse(S.be)
    S.be = full(S.be);
end

% equality constraints are treated as double-sided inequalities
G = [S.A; S.Ae];
g2 = [S.b; S.be];
% lower bound for the left hand side of inequalities g1 <= A*x<= g2 is
% set as very low number to achieve dimension match
if ~S.test
    g1 = [-MPTOPTIONS.infbound*ones(S.m,1); S.be];
else
    g1 = [-1e9*ones(S.m,1); S.be];
end
%The integer argument nWSR specifies the maximum number of working set
%recalculations to be performed during the initial homotopy (on output it contains the number
%of working set recalculations actually performed!)
nWSR = 5*(S.m+S.n); % this value is suggested by qpOASES manual


% call qpOASES
[R.obj, R.xopt, lambda, status]=qpOASES(S.H, S.f, G, S.lb,...
    S.ub, g1, g2, nWSR, S.x0);

% extract multipliers
if ~isempty(S.lb)
    if S.test
        activelb = (R.xopt < S.lb + 1e-4 );
    else
        activelb = (R.xopt < S.lb + MPTOPTIONS.rel_tol );
    end
else
    activelb = false(S.n,1);
end
if ~isempty(S.ub)
    if S.test
        activeub = (R.xopt > S.ub - 1e-4 );
    else
        activeub = (R.xopt > S.ub - MPTOPTIONS.rel_tol );
    end
else
    activeub = false(S.n,1);
end
R.lambda.lower = lambda(1:S.n);
R.lambda.lower(~activelb) = 0;
R.lambda.upper = -lambda(1:S.n);
R.lambda.upper(~activeub) = 0;
R.lambda.ineqlin = -lambda(S.n+1:S.n+S.m);
R.lambda.eqlin = -lambda(S.n+S.m+1:S.n+S.m+S.me);


switch status
    case 0
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case 1
        R.how = ['Maximum number of working sets ', num2str(nWSR),' reached.'];
            %'You may change this value in S.options.nWSR.'];
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end

    otherwise
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end

end
