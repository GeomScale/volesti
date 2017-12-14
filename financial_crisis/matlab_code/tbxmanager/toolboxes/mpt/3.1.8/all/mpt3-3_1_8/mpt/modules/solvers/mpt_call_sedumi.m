function R = mpt_call_sedumi(S)
%
%  MPT_CALL_SEDUMI: A gateway function to SEDUMI solver (without errorchecks) 
%  ===========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_sedumi(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements calls to solve LP directly and QP via transformation
%  to second order cone. For LP in a form 
%                                          T                  
%                                   min   f x             (1) 
%                                                             
%                                  s.t.   Ax <= b         (2) 
%                                                             
%                                         A x = b         (3) 
%                                          e     e            
%     we need to get a following form accepted by SEDUMI 
%                                          T                  
%                                  min   f  x             (4) 
%                                         f                   
%                                                             
%                                 s.t.   A x >= b         (5) 
%                                         f      f            
%                                        x >= 0           (6) 
%     This can be achieved by introducing variables x^+ >= 0, x^- >= 0, and y >= 0 
%                                         +    -  
%                                    x = x  - x   
%                                                 
%                                    y = -Ax+b    
%                                                 
%     and putting them in one vector z = [x^+,  x^-,  y]. The LP to be solved by
%  SEDUMI is formed by 
%                                       (         )            
%                                 f  =  ( f -f 0  )        (7) 
%                                  f    (         )            
%                                       ( -A  A  -I  )         
%                                       (            )         
%                                 A  =  ( A  -A  0   )     (8) 
%                                  f    (  e   e     )         
%                                       ( -b  )                
%                                       (     )                
%                                 b  =  ( b   )            (9) 
%                                  f    (  e  )                
%     and solved in z  variables. For QP of the form 
%                                        1 T    T              
%                                  min   -x Hx f x        (10) 
%                                        2                     
%                                 s.t.   Ax <= b          (11) 
%                                                              
%                                        A x = b          (12) 
%                                         e     e              
%     we need to impose constraints a in quadratic cone K  and to express the above
%  problem as 
%                                          T                  
%                                   min   c x            (13) 
%                                          n                  
%                                                             
%                                  s.t.   A x = b        (14) 
%                                          n     n            
%                                         x in K         (15) 
%     Since QP is convex, we can express (??)-(??)  in an epigraph form 
%                                                                 
%                              min   t                       (16) 
%                                                                 
%                             s.t.   Ax <= b                 (17) 
%                                                                 
%                                    A x = b                 (18) 
%                                     e     e                     
%                                    1 T      T                   
%                                    -x Hx + f x <= t        (19) 
%                                    2                            
%                                                            (20) 
%     over quadratic constraints (??). From the literature for convex programming
%  follows that the quadratic constraints (??)  can be written as conic constraints
%   {                                        (             ||   Qx     ||  ) }     
%                                           
%   {                                        (   T         ||    T     ||  ) }     
%                                           
%   {           ( T T      T            )    (1-f x+t      || 1+f x-t  ||  ) }     
%                                           
%    { (x,t)  |  (x Q Qx + f x -t  <=  0 ) =  (         <=  ||          ||  ) }   
%                                        (21)
%   {           (                       )    (-------      || -------  ||2 ) }     
%                                           
%   {                                        (   2         ||    2     ||  ) }     
%                                           
%     where the matrix Q  is a Cholesky factor of 0.5H = Q^TQ. This equivalence
%  allows us to rewrite the epigraph form of QP (??)-(??)to a primal form (??)-(??)
%   accepted by SEDUMI. Equality and inequality constraints are treated the same
%  way as in LP case, i.e. by introducing the new variables x = x^+ - x^-, x^+ >=
%  0, x^- >= 0, and y=-Ax+b, y >= 0. Moreover, to express conic constraints we need
%  two more variables v in R^1, u in R^n+1 
%                                         T         
%                                      1-f x+t      
%                                  v =              
%                                      -------      
%                                         2         
%                                      (   Qx     ) 
%                                      (    T     ) 
%                                      ( 1+f x-t  ) 
%                                  u = (          ) 
%                                      ( -------  ) 
%                                      (    2     ) 
%     Collecting all variables to one column vector z = [t,  x^+,  x^-,  y,  v,  u]
%   the linear equality constraints (??)  can be expressed in z  variable as  A_nz
%  = c_c  where 
%                                ( 0  -A   A  -I 0  0   ) 
%                                (                      ) 
%                                ( 0  A   -A  0  0  0   ) 
%                                (     e    e           ) 
%                                ( 0   Q  -Q  0  0  -I  ) 
%                                (     T    T           ) 
%                                ( -1 f   -f            ) 
%                           A  = ( --         0  0  -1  ) 
%                            n   ( 2  --  ---           ) 
%                                (    2    2            ) 
%                                (      T  T            ) 
%                                ( 1  -f  f             ) 
%                                ( -          0  -1 0   ) 
%                                ( 2  --- --            ) 
%                                (     2  2             ) 
%                                ( -b  )                  
%                                (     )                  
%                                ( b   )                  
%                                (  e  )                  
%                                ( 0   )                  
%                           b  = ( -1  )                  
%                            n   ( --  )                  
%                                ( 2   )                  
%                                ( -1  )                  
%                                ( --  )                  
%                                ( 2   )                  
%     The objective function in (??)  is composed as c_n^T = [1,  0,  0,  0,  0, 
%  0]. The individual constraints K  in (??)  are given in z  vector as follows: t 
%  is free variable, x^+, x^-, y  are restricted to be nonnegative and u, v  form
%  the conic constraint v>=||u||_2.
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
%                         solved (only LP and QP are allowed).     
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
%          R.lambda   Lagrangian multipliers                   
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

if ~any(strcmpi(S.problem_type,{'LP','QP'}))
    error('mpt_call_clp: SEDUMI solver does not solve %s problems!',S.problem_type);
end

% merge inequality constraints
A=S.A;
b=S.b;

% detect Inf boundaries
if S.test
    ilb = (S.lb==-Inf) | (S.lb<=-1e6);
    iub = (S.ub==Inf)  | (S.ub>=1e6);
else
    ilb = (S.lb==-Inf) | (S.lb<=-MPTOPTIONS.infbound);
    iub = (S.ub==Inf)  | (S.ub>=MPTOPTIONS.infbound);
end
% store kept rows
kept_rows.lb = find(~ilb);
kept_rows.ub = find(~iub);
if any(~ilb)
    % put ones at the positions where there is lb/ub
    Alb = zeros(nnz(~ilb),S.n);
    Alb(:,~ilb) = -eye(nnz(~ilb));
    A = [A; Alb];
    b = [b; -S.lb(~ilb)];
end
if any(~iub)
    Aub = zeros(nnz(~iub),S.n);
    Aub(:,~iub) = eye(nnz(~iub));
    A = [A; Aub];
    b = [b; S.ub(~iub)];
end
m = size(A,1);

if strcmpi(S.problem_type,'LP')
   
    % convert LP to form accepted by Sedumi
    %     min   c'*x
    % s.t.:   A*x = b
    %          x >= 0
   
    % Equality constraints Ae*x = be are extended in variables xp, xm
    % x = xp - xm
    % Inequality constraies A*x <= b are replaced by slack variable y to be
    % positive, -A*x-y = -b, y>=0
    
    % the final form from equalities
    % [Ae -Ae 0]*[xp; xm; y] = be
   
    % the final form from inequalities
    % [-A A -I]*[xp; xm; y] = -b
    %   xp, xm, y >= 0
       
    An = [-A A -eye(m)];
    bn = -b;
    if S.me>0
        An = [S.Ae -S.Ae zeros(S.me,m); An];
        bn = [S.be; bn];
    end
    
    % set dimension of nonnegativity constraints
    K.l = S.n+S.n+m;
    
    % objective function
    cn = [S.f; -S.f; zeros(m,1)];
    
    % call sedumi    
    if S.test
        % use default settings, with fid=0 (no verbosity)
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,struct('fid',0));
    else
        % use global MPT settings
        opts = MPTOPTIONS.modules.solvers.sedumi;
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,opts);
    end

    % recover original variables
    if ~isempty(z)
       R.xopt = full(z(1:S.n)-z(S.n+1:2*S.n));
    else
       R.xopt = zeros(S.n,1); 
    end
            
else
    % for QP we need to create second order cone
    %
    %     min t
    % s.t.   Ae*x = be
    %        A*x <= b
    %   x'*(H/2)*x + f'*x <= t
    %
    % Linear inequalities and equalities are rewritten in xp, xm, y
    % variables as in LP.
    % The quadratic constraint can be written as
    %
    %  ||      Qx            ||
    %  || 0.5*(1 + f'*x -t)  ||_2   <= 0.5*(1-f'*x+t)
    %  
    % which is passed to Sedumi via substitution
    %  u = [Q*x; 0.5*(1+f'*x-t)];
    %  v = [0.5*(1-f'*x+t)];
    % 
    % where H/2 = Q'*Q
    
    % factor H/2 = Q'*Q
    if S.test
        Q = cholinc(sparse(S.H/2),1e-8);
    else
        Q = cholinc(sparse(S.H/2),MPTOPTIONS.abs_tol);
    end
    % objective function in variables z=[t(1); xp(n); xm(n); y(m); v(1); u(n+1)]
    cn =[1; zeros(3*S.n+m+2,1)];
    
    % inequalities and equalities are extended in z variables as in LP case
    An = [zeros(m,1) -A A -eye(m) zeros(m,S.n+1) zeros(m,1)];
    bn = -b;
    if S.me>0
        An = [zeros(S.me,1) S.Ae -S.Ae zeros(S.me,m) zeros(S.me,S.n+1) zeros(S.me,1); An];
        bn = [S.be; bn];
    end
    
    % add second order cone constraints
    An = [An;
        zeros(S.n,1) Q -Q zeros(S.n,m) zeros(S.n,1) -eye(S.n,S.n+1);
       -0.5   S.f'/2  -S.f'/2  zeros(1,m) 0 [zeros(1,S.n) -1];
        0.5  -S.f'/2   S.f'/2  zeros(1,m) -1 zeros(1,S.n+1)];
    bn = [bn; zeros(S.n,1); -0.5; -0.5];
    
    % t, xp, xm, y >= 0
    % u, v belong to quadratic cone
    K.f = 1; % t is a free variable
    K.l = S.n+S.n+m; % nonnegative variables
    K.q = S.n+1+1; % v, u belong to quadratic cone v >= norm(u)
    
    % call sedumi in sparse format
    if S.test
        % use default settings, with fid=0 (no verbosity)
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,struct('fid',0));
    else
        % use global MPT settings
        opts = MPTOPTIONS.modules.solvers.sedumi;
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,opts);
    end
        
    % recover original variables in full format
    if ~isempty(z)
        R.xopt = full(z(2:S.n+1)-z(S.n+2:2*S.n+1));
    else
        R.xopt = zeros(S.n,1);
    end
        
end

% recover multipliers for both LP/QP case because the order of constraints
% remains the same (equalities, inequalities, lb, ub)
R.lambda.ineqlin = full(dz(S.me+1:S.me+S.m));
R.lambda.eqlin = -full(dz(1:S.me));
if ~isempty(S.lb)
    R.lambda.lower = zeros(S.n,1);
    R.lambda.lower(kept_rows.lb) = full(dz(S.me+S.m+1:S.me+S.m+numel(kept_rows.lb)));
else
    R.lambda.lower = zeros(S.n,1);
end
if ~isempty(S.ub) && isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = full(dz(S.me+S.m+1:S.me+S.m+numel(kept_rows.ub)));
elseif ~isempty(S.ub) && ~isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = full(dz(S.me+S.m+numel(kept_rows.lb)+1:S.me+S.m+numel(kept_rows.lb)+numel(kept_rows.ub)));
else
    R.lambda.upper = zeros(S.n,1);
end

% if the field is missing, we assume infeasibility
if ~isfield(status,'numerr')
    status.numerr = 3;
end
if ~isfield(status,'pinf')
    status.pinf = 0;
end
if ~isfield(status,'dinf')
    status.dinf = 0;
end


% check numerr first
switch status.numerr
    case 0
        R.how = 'ok';      
    case 1
        R.how = 'ok, but with numerical problems';
    case 2
        R.how = 'numerical problems';
    otherwise
        R.how = 'unknown problem';
end

% check primal, dual feasiblity
if status.pinf==status.dinf && any(status.numerr==[0,1])
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
elseif status.pinf==1,
    R.how = [R.how,', primal infeasible'];
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
elseif status.dinf==1
    R.how = [R.how,', dual infeasible'];
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
else
    if S.test
        R.exitflag = -1;
    else
        R.exitflag = MPTOPTIONS.ERROR;
    end
end

R.obj = S.f'*R.xopt;
if ~isempty(S.H)
    R.obj = R.obj + 0.5*R.xopt'*S.H*R.xopt;
end
