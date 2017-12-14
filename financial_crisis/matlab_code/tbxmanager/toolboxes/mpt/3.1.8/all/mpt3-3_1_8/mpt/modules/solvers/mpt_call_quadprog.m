function R = mpt_call_quadprog(S)
%
%  MPT_CALL_QUADPROG: A gateway function to QUADPROG solver (without errorchecks) 
%  ===============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_quadprog(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements call to QUADPROG solver based on formulation from Opt
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
%     which is passed to QUADPROG solver directly. For LP QUADPROG solver passes
%  the data to LINPROG automatically.
%  
%  INPUT
%  -----
%     
%        
%          S              Structure of the Opt class               
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

if ~any(strcmpi(S.problem_type,{'QP'}))
    error('mpt_call_quadprog: QUADPROG solver does not solve %s problems!',S.problem_type);
end

% overwrite default settings
if S.test
    try
        % MOSEK
        % I really hate doing this. Hey, MOSEK, stop messing with linprog!
        options = mskoptimset('quadprog');
        options.Display = 'off';
    catch
        options=optimset(optimset('quadprog'),'Display','off');
    end
else
    options=MPTOPTIONS.modules.solvers.quadprog;
end

% make sure the Hessian is symmetric.options=optimset(optimset('quadprog'),'Display','off','LargeScale','off',...
if norm(S.H-S.H', Inf) < 1e-10,
    % but we only remove numerical noise if the hessian is only
    % "slightly" wrong. for clearly non-symmetrical hessians we still
    % let quadprog to display a proper warning
    S.H = (S.H + S.H')*0.5;
end

% direct call to quadprog
[R.xopt,R.obj,exitflag,OUTPUT,R.lambda]=quadprog(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,S.x0,options);
if exitflag>0 %then QUADPROG converged with a solution X.
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
else
    % ==0 then the maximum number of iterations was exceeded (only occurs
    %     with large-scale method).
    % < 0 then the problem is unbounded, infeasible, or
    %     QUADPROG failed to converge with a solution X.
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
end
%R.lambda=[lambdav.ineqlin; lambdav.eqlin];

end
