function R = mpt_call_qpspline(S)
%
%  MPT_CALL_QPSPLINE: A gateway function to QPspline solver (without errorchecks) 
%  ===============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_qpspline(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements call to QPspline solver based on formulation from Opt
%  class. Only QP problems are supported with positive definite Hessian. It is
%  assumed that QP entering this function is of the form 
%                                        1  T    T             
%                                  min   - x Hx+f x        (1) 
%                                        2                     
%                                 s.t.   Ax <= b           (2) 
%                                                              
%                                        A x = b           (3) 
%                                         e     e              
%     which must be transformed to 
%                                      1  T      T              
%                                min   - x H x+f  x         (4) 
%                                      2    f   f               
%                                                               
%                               s.t.   l <= A x <= u        (5) 
%                                            f                  
%                                                           (6) 
%     which accepts QPspline. The lower bound l  is always set as
%  -MPTOPTIONS.infbound. If QP contains equality constraints, these are removed
%  first. It is required that the system of linear equations A_ex=b_e  is
%  consistent, i.e. no linearly dependent rows are found and the number of
%  equalities is strictly less than number of variables. The principle is based on
%  factorizing equality constraints A_ex=b_e  in basic  x_Bc and non-basic
%  variables x_Nc, i.e. 
%                                    (              )
%                               A  = ( A     A      )
%                                e   (  e,Bc  e,Nc  )
%    which gives 
%                                                      
%                             A    x   + A    x   = b  
%                              e,Bc Bc    e,Nc Nc    e 
%     where the index sets Bc, Nc denote the columns from which factored system is
%  built. The factored submatrix A_e,Bc  must be invertible in order to express
%  basic variables as a function of non-basic variables, i.e. 
%                                   -1                -1      
%                       x   = -A      A    x   + A      b     
%                        Bc     e,Bc   e,Nc Nc    e,Bc   e,Bc 
%     With the substitution 
%                                           -1      
%                                 C = -A      A     
%                                       e,Bc   e,Nc 
%     and 
%                                          -1      
%                                 D = A      b     
%                                      e,Bc   e,Bc 
%    the relation between basic and non-basic variables is simplified to 
%                                                      
%                               x   = Cx   + D      (7)
%                                Bc     Nc             
%     The above QP problem (??)-(??)  can be expressed only in non-basic variables
%  x_Nc  as follows: 
%                                  1    T      T                    
%                            min   - x   Hx  +f x   + c         (8) 
%                                  2  Nc   Nc    Nc                 
%                                                                   
%                           s.t.   Ax   <= b                    (9) 
%                                    Nc                             
%                                                              (10) 
%     where 
%                          T           T                                    
%                   H  =  C H     C + C H      + H     C + H           (11) 
%                            Bc,Bc       Bc,Nc    Nc,Bc     Nc,Nc           
%                          T           T            T       T               
%                   f  =  D H     C + D H      + f   C + f             (12) 
%                            Bc,Bc       Nc,Bc    Bc      Nc                
%                         1 T             T                                 
%                   c  =  -D H     D + f   D                           (13) 
%                         2   Bc,Bc     Bc                                  
%                                                                           
%                   A  =  A  C+A                                       (14) 
%                          Bc   Nc                                          
%                                                                           
%                   b  =  b - A  D                                     (15) 
%                              Bc                                           
%     Original solution to QP problem (??)-(??)  can be obtained via relation (??).
%  
%  
%  INPUT
%  -----
%     
%        
%          S              structure of the Opt class               
%                         Class: struct                            
%          S.H            Quadratic part of the objective function 
%                         which is strictly convex Hsucc 0.        
%                         Class: double                            
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
%          S.lb           Lower bound for the variables x >= l_b.  
%                         Class: double                            
%                         Default: []                              
%          S.ub           Upper bound for the variables x <= u_b.  
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

if ~strcmpi(S.problem_type,'QP')
    error('mpt_call_clp: QPspline solver does not solve %s problems!',S.problem_type);
end

% substitute
H = S.H;
f = S.f;
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

% store merged constraints
Am = A;
bm = b;


% remove equality constraints if any
Ae = S.Ae;
be = S.be;
me = S.me;
kept_rows.eq = 1:S.me;
if me>0
    % factorize Ae to get
    %  Ae(Br,Bc)*x(Bc) + Ae(Br,Nc)*x(Nc) = be(Br) % must be invertible mapping
    
    % check rank of Ae
    if S.test
        re = rank(Ae);
    else
        re = rank(Ae,MPTOPTIONS.rel_tol);
    end
    % check consistency
    if S.test
        rc = rank([Ae be]);
    else
        rc = rank([Ae be],MPTOPTIONS.rel_tol);
    end
    
    % for underdetermined system check linearly dependent rows
    if re<me || rc<me
        
        % if the right hand side is not linearly dependent, infeasible solution
        if rc>re
            R.xopt = zeros(S.n,1);
            R.obj = 0;
            R.lambda.ineqlin = [];
            R.lambda.eqlin = [];
            R.lambda.lower = [];
            R.lambda.upper = [];
            if S.test
                R.exitflag = 2;
            else
                R.exitflag = MPTOPTIONS.INFEASIBLE;
            end
            R.how = 'infeasible';
            return
        end
        
        while me ~= re
            % find linearly dependent rows
            [~,~,p] = lu(sparse([Ae be]),'vector');
            rd = p(re+1:end);
            
            % remove linearly dependent rows
            Ae(rd,:) = [];
            be(rd) = [];
            me = me-length(rd);
            kept_rows.eq(rd) = [];
            
            if S.test
                re = rank(full([Ae be]));
            else
                re = rank(full([Ae be]),MPTOPTIONS.rel_tol);
            end
        end
        
    end
    
    [Le,Ue,pe,qe] = lu(sparse(Ae),'vector');
    Br = pe(1:re); Bc = qe(1:re);
    Nr = pe(re+1:end); Nc = qe(re+1:end);
    % Nr must be empty -> otherwise there are depentent rows in S.Ae
    % which must be removed
    
    % substitute x(Bc) = C*x(Nc) + D
    Aebn = S.Ae(Br,Nc);
    %iAebb = inv(S.Ae(Br,Bc));
    beb = S.be(Br);
    
    % use factorized solution to compute C
    % C = -S.Ae(Br,Bc)\Aebn;
    Cl = -linsolve(full(Le),Aebn,struct('LT',true));
    C = linsolve(full(Ue(:,1:re)),Cl,struct('UT',true));
    
    % use factorized solution to compute D
    % D = S.Ae(Br,Bc)\beb;
    Dl = linsolve(full(Le),beb,struct('LT',true));
    D = linsolve(full(Ue(:,1:re)),Dl,struct('UT',true));
    
    Abc = A(:,Bc); Anc = A(:,Nc);
    
    % modify inequality constraints
    %A = -Abc*iAebb*Aebn + Anc;
    A = Abc*C + Anc;
    %b = b - Abc*iAebb*beb;
    b = b - Abc*D;
    
    % modify cost
    %H = S.H(Nc,Nc) + Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*Aebn - Aebn'*iAebb'*S.H(Bc,Nc) - S.H(Nc,Bc)*iAebb*Aebn;
    H = S.H(Nc,Nc) + C'*S.H(Bc,Bc)*C + C'*S.H(Bc,Nc) + S.H(Nc,Bc)*C;
    %f = S.H(Nc,Bc)*iAebb*beb - Aebn'*iAebb'*S.f(Bc) + S.f(Nc) - Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*beb;
    f = S.H(Nc,Bc)*D + C'*S.f(Bc) + S.f(Nc) + C'*S.H(Bc,Bc)*D;
    
end

% actual dimensions
[m,n] = size(A);

% transformation to form
% min 0.5 x'Mx - b'x
% st:  l< Ax <u

% artificial lower bound
if ~S.test
    l = -MPTOPTIONS.infbound*ones(m,1);
else
    l = -1e9*ones(m,1);
end

if S.test
    [R.xopt, lambda, R.obj, exitflag, R.how, iter, time] =  QPspline(H, f, l, A, b);
else
    [R.xopt, lambda, R.obj, exitflag, R.how, iter, time] =  QPspline(H, f, l, A, b, MPTOPTIONS.modules.solvers.qpspline);
end


% extract Lagrange multipliers
R.lambda.ineqlin = -lambda(1:S.m);
if ~isempty(S.lb)
    R.lambda.lower = zeros(S.n,1);
    R.lambda.lower(kept_rows.lb) = -lambda(S.m+1:S.m+numel(kept_rows.lb));
else
    R.lambda.lower = zeros(S.n,1);
end
if ~isempty(S.ub) && isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = -lambda(S.m+1:S.m+numel(kept_rows.ub));
elseif ~isempty(S.ub) && ~isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = -lambda(S.m+numel(kept_rows.lb)+1:S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
else
    R.lambda.upper = zeros(S.n,1);
end


% if there were equalities, map back to original solution
if S.me>0
    xopt = zeros(S.n,1);
    xopt(Nc) = R.xopt;
    xopt(Bc) = C*R.xopt + D;
    R.xopt = xopt;
    % solve overdetermined system to get multipliers for equalities
    % H*x + f + Am'*lambda_ineq + Ae'*lambda_eq = 0
    lambda_eq = zeros(S.me,1);
    lambda_eq(kept_rows.eq) = -Ae'\(S.H*R.xopt + S.f - Am'*lambda);
    
    % extend multipliers
    R.lambda.eqlin = lambda_eq;
    R.obj = 0.5*R.xopt'*S.H*R.xopt + S.f'*R.xopt;
else
    R.lambda.eqlin = []; 
end


% correct multipliers for lower/upper bounds
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
R.lambda.lower(~activelb) = 0;
R.lambda.upper(~activeub) = 0;

% recast exitflag to MPT type
switch exitflag
    case 1
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case -1
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end

    case {-2,-3,-4}
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end        
end
    
    

end
