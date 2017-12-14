function R = mpt_call_qpc(S)
%
%  MPT_CALL_QPC: A gateway function to QPC solver (without errorchecks) 
%  =====================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_qpc(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements call to QPC solver based on formulation from Opt
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
%     which is passed to QPC solver directly. Sparse inputs are converted to full
%  if needed. QPC offers two types of algorithms to solve QP/LP. For an interior
%  point method specify in the field "solver" a string "qpip". Otherwise active set
%  method is chosen by default.
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
%          S.solver       Specific routine to be called of QPC. To 
%                         call interior point method, specify      
%                         "qpip". To call active set method,       
%                         specify "qpas" or leave empty.           
%                         Class: char                              
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
%                     appeared through optimization            
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

if ~any(strcmpi(S.problem_type,{'LP','QP'}))
    error('mpt_call_qpc: QPC solver does not solve %s problems!',S.problem_type);
end

% for LP set H=0
if isempty(S.H)
    S.H = zeros(S.n);
end

% QPC does not accept sparse matrices
if issparse(S.H)
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

show_display = 0;

if any(strcmpi(S.solver,{'qpip','qpc.qpip'}))

    if S.test        
        % primal-dual predictor-corrector algorithm    
        [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display);
    else
        mu = MPTOPTIONS.modules.solvers.qpip.mu;
        method = MPTOPTIONS.modules.solvers.qpip.method;
        [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display,mu,method);
    end
else 
    % dual active-set algorithm (default)
    [x, ef, L] = qpas(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display);
end

% test for NaN, or Inf
if any(isnan(x)) || any(isinf(x))
    x = zeros(S.n,1);
    ef = -1;
end

% must check feasibility
if ef==0 
    if S.m>1
        if S.test
            ve = any(S.A*x>S.b+1e-5);
        else
            ve = any(S.A*x>S.b+MPTOPTIONS.rel_tol);
        end
    else
        ve = false;
    end
    if S.me>1
        if S.test
            vq = ( norm(S.Ae*x-S.be,Inf) > 1e-5 );
        else
            vq = ( norm(S.Ae*x-S.be,Inf) > MPTOPTIONS.rel_tol );
        end
    else
        vq = false;
    end
    if any(ve) || vq
        % retry with interior point method
        if S.test
            [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display);
        else
            mu = MPTOPTIONS.modules.solvers.qpip.mu;
            method = MPTOPTIONS.modules.solvers.qpip.method;
            [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display,mu,method);
        end
        % check again feasibility
        if ef==0
            if S.m>1
                if S.test
                    ve = any(S.A*x>S.b+1e-5);
                else
                    ve = any(S.A*x>S.b+MPTOPTIONS.rel_tol);
                end
            else
                ve = false;
            end
            if S.me>1
                if S.test
                    vq = ( norm(S.Ae*x-S.be,Inf) > 1e-5 );
                else
                    vq = ( norm(S.Ae*x-S.be,Inf) > MPTOPTIONS.rel_tol );
                end
            else
                vq = false;
            end
            if any(ve) || vq                
                % infeasible
                ef = -1;
            end
        end
    end
end

R.xopt = x;
R.obj = 0.5*x'*S.H*x + S.f'*x;
R.lambda.ineqlin = L.inequality;
R.lambda.eqlin = L.equality;
R.lambda.lower = L.lowerbound;
R.lambda.upper = L.upperbound;

if ef==0
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
elseif ef==-1
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end    
else
    R.how = 'unknown (possibly infeasible)';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
    % since we don't know what the other flags mean, throw an error
    %error('mpt_call_qpc: other error code in "exitflag" from qpc.')
end
