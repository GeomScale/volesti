function R=mpt_solve(Sin, param_value)
%
%  MPT_SOLVE: A gateway function to solve non-parametric optimization problems
%  ===========================================================================
%  (without errorchecks) 
%  ======================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_solve(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function is the main routine for fast calls for solving non-parametric
%  optimization problems in MPT. In fact, it is a subroutine of Opt as a part of
%  solve method. The Opt class serves as general wrapper for preprocessing the data
%  involved in optimization, including necessary error checks. Once the data are
%  valid, then are passed to mpt_solve function that calls the appropriate solver.
%  It is assumed that QP/LP/MIQP/MILP and entering this function (for LP/MILP H=0)
%  is of the form 
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
%    which is given by strings in vartype field. The matrices H, A, A_e, and
%  vectors f, b, b_e, lb, ub are the problem data, then x  are the decision
%  variables. The LCP must be given as: 
%                                    w - Mz= q         (6) 
%                                       w >= 0         (7) 
%                                       z >= 0         (8) 
%                                      T                   
%                                     w z = 0          (9) 
%                                                          
%     where the matrices M, Q, and vectors q  are the problem data, z, w  are the
%  decision variables to be determined. The function mpt_solve processes the
%  problem data and passes it to the appropriate solver. The particular solver can
%  be specified by providing solver name or it is selected automatically from the
%  list of available solvers. Based on the solver name, the appropriated function
%  is called: 
%    
%     - CDD solver is installed by default, and can be called via mpt_call_cdd.
%     Solves LP problems. 
%     - CLP solver is installed by default, and can be called via mpt_call_clp.
%     Solves LP/QP problems. 
%     - CPLEX is a commercial solver and must be installed additionally. It can be
%     called via mpt_call_cplex.  CPLEX solves LP/QP/MILP/MIQP problems. 
%     - GLPK solver is installed by default, and can be called via mpt_call_glpk.
%     It solves LP/QP/MILP problems. 
%     - GUROBI is a commercial solver and must be installed additionally. It can be
%     called via mpt_call_gurobi.  It solves LP/QP/MILP/MIQP problems. 
%     - LCP is a default solver, and can be called via mpt_call_lcp. It can be used
%     to solve LP/QP and LCP problems. 
%     - LINPROG is Matlab LP solver, and can be called via mpt_call_linprog. 
%     - QUADPROG is Matlab QP solver, and can be called via mpt_call_quadprog. 
%     - NAG is a commercial solver and must be installed additionally. It can be
%     called via mpt_call_nag.  It solves LP/QP problems. 
%     - QPC is LP/QP solver that need to be installed additionally. It can be
%     called via mpt_call_qpc. 
%     - QPOASES is LP/QP solver that is installed by default. It can be called via
%     mpt_call_qpoases. 
%     - QPSPLINE is a QP solver for strictly convex problems and can be called via
%     mpt_call_qpspline. 
%     - SEDUMI is a semidefinite solver for general convex problems and can be
%     called via mpt_call_sedumi. 
%  
%  
%  INPUT
%  -----
%     
%        
%          S              structure of the Opt class               
%                         Class: struct                            
%          S.H            quadratic part of the objective function 
%                                                                  
%                         Class: double                            
%                         Default: []                              
%          S.f            linear part of the objective function    
%                         Class: double                            
%          S.A            linear part of the inequality            
%                         constraints Ax <= b                      
%                         Class: double                            
%          S.b            right hand side of the inequality        
%                         constraints Ax <= b                      
%                         Class: double                            
%          S.Ae           linear part of the equality constraints  
%                         A_ex=b_e                                 
%                         Class: double                            
%                         Default: []                              
%          S.be           right hand side of the equality          
%                         constraints A_ex=b_e                     
%                         Class: double                            
%                         Default: []                              
%          S.lb           lower bound for the variables x >= lb    
%                         Class: double                            
%                         Default: []                              
%          S.ub           upper bound for the variables x <= ub    
%                         Class: double                            
%                         Default: []                              
%          S.M            Positive semi-definite matrix defining   
%                         LCP.                                     
%                         Class: double                            
%                         Default: []                              
%          S.q            Right hand side vector defining LCP.     
%                         Class: double                            
%                         Default: []                              
%          S.n            problem dimension (number of variables)  
%                         Class: double                            
%          S.m            number of inequalities in Ax <= b        
%                         Class: double                            
%          S.me           number of equalities in A_ex=b_e         
%                         Class: double                            
%          S.problem_type a string specifying the problem to be    
%                         solved                                   
%                         Class: char                              
%          S.vartype      A string specifying the type of          
%                         variable. Supported characters are C     
%                         (continuous), I (integer), B (binary), N 
%                         (semi-integer), S (semi-continuous).     
%                         Example: First variable from three is    
%                         binary, the rest is continuous:          
%                         S.vartype='BCC';                         
%                         Class: char                              
%          S.solver       S string specifying which solver should  
%                         be called.                               
%                         Class: char                              
%          S.test         Call (false) or not to call (true) MPT   
%                         global settings                          
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
%          R.xopt     Optimal solution                         
%                     Class: double                            
%          R.obj      Objective value                          
%                     Class: double                            
%          R.lambda   Lagrangian multipliers                   
%                     Class: double                            
%          R.exitflag An integer value that informs if the     
%                     result was feasible (1), or otherwise    
%                     (different from 1)                       
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
%     Opt
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

% speed-hack
if isstruct(Sin) && isfield(Sin, 'quicklp') && Sin.quicklp
    S = Sin;
    [S.m, S.n] = size(Sin.A);
    S.me = size(Sin.Ae, 1);
    S.f = S.f(:);
    S.solver = MPTOPTIONS.lpsolver;
    S.problem_type = 'LP';
    S.x0 = zeros(S.n,1);
    S.test = false;
    S.d = [];
    S.c = 0;
    S = sub_fix_nan(S);
    R = sub_callsolver(S, MPTOPTIONS);
    return

elseif isstruct(Sin) && isfield(Sin, 'quickqp') && Sin.quickqp
    S = Sin;
    [S.m, S.n] = size(Sin.A);
    S.me = size(Sin.Ae, 1);
    S.f = S.f(:);
    S.solver = MPTOPTIONS.qpsolver;
    S.problem_type = 'QP';
    S.x0 = zeros(S.n,1);
    S.test = false;
    S.d = [];
    S.c = 0;
    S = sub_fix_nan(S);
    R = sub_callsolver(S, MPTOPTIONS);
    return
end

narginchk(1, 2);

% assuming that the problem is defined in a structure S of Opt class
% 
% LP/QP/pLP/pQP variables
% J(th) = min 0.5*x'*H*x + (pF*th+f)'*x + th'*Y*th + C*th + c
%         s.t.  A*x <= b  + pB*th
%               Ae*x = be + pE*th
%               lb  <= x <= ub
%               Ath*th <= bth
%
% LCP variables
% w - M*z = q + Q*th, w,z  >= 0, w'*z = 0
%
%
% S should contain fields:
%   
%   H - quadratic part of the objective function
%   f - linear part of the objective function
%   pf - parametric linear part of the objective function
%   A - constraints inequality set in A*x <= b
%   b - rhs of inequalities
%   pB - parametric rhs of inequalities
%   Ae - equality constraint set in Ae*x = be
%   be - rhs of equality constraints
%   pE - parametric rhs of equality constraints
%   lb - lower bound on optimization variables x
%   ub - upper bound on optimization variables x
%   Ath - inequality set for parametric variables th in Ath*th <= bth
%   bth - rhs of the inequality set for parameters th
%   Y - quadratic term for parameters in the cost function
%   C - linear term for parameters in the cost function
%   c - constant term in the objective function
%   M - LCP
%   q - LCP
%   Q - LCP
%   n - Problem dimension
%   m - Number of inequalities
%   me - Number of equality constraints
%   d - Number of parameters
%   solver - solver to be used

arg_set = {'H','f','pF','A','b','pB','Ae','be','pE','lb','ub','Ath','bth', ...
    'Y','C','c','M','q','Q','n','m','me','d','solver'};

% added fields: 
%   x0 - initial conditions
%   vartype - specifies which variables are binaries and which continuous
%   problem_type - 'LP', 'QP', ...
%   test - solve a trivial problem to test if the solver works - used in
%         the initialization phase by "mpt_detect_solvers"
%   routine - routine to choose for LCP solver

opt_arg_set = {'x0','vartype','problem_type','test','routine'};

% construct new object with empty fields
S = sub_initialize_args(arg_set, opt_arg_set);

if isa(Sin,'Opt')
    % all fields are present and checked
    % create new structure by copying all fields
    fnames = fieldnames(Sin);
    for i=1:numel(fnames)
        S.(fnames{i}) = Sin.(fnames{i});
    end
        
else    
    % assign values from Sin
    fs = fieldnames(Sin);
    for i=1:length(fs)
        S.(fs{i}) = Sin.(fs{i});
    end
    
    % this field is used at the initialization of the toolbox, to actually test
    % if the solver works
    % used in "mpt_detect_solvers" routine
    % not used in parametric case, because here we do not test solvers
    if isempty(S.test)
        S.test = false;
    end
    
    % global options
    if isempty(MPTOPTIONS) && ~S.test
        MPTOPTIONS = mptopt;
    end
    
    % checking for required arguments
    if isempty(S.M) % LP,QP       
        if isempty(S.A) && isempty(S.f) && isempty(S.Ae)
            error('mpt_solve: One of the fields "A", "Ae", or "f" must be nonempty!');
        end
    else % LCP
        if isempty(S.q)
            error('mpt_solve: Field "q" must not be empty.')
        end
        S.n = length(S.q);
        S.m = 0;
        S.d = size(S.Q,2);
    end
    
    
    % some arguments must not be empty
    if isempty(S.n)
        if ~isempty(S.A)
            nA = size(S.A,2);
        else
            nA = 0;
        end
        if ~isempty(S.Ae)
            nAe = size(S.Ae,2);
        else
            nAe = 0;
        end
        if ~isempty(S.f)
            nf = max(size(S.f));
        else
            nf = 0;
        end
        S.n = max([nA nAe nf]);
    end
    if isempty(S.m)
        S.m = size(S.A,1);
    end
    if isempty(S.me)
        if ~isempty(S.Ae)
            S.me = size(S.Ae,1);
        else
            S.me = 0;
        end
    end
    if isempty(S.d)
        if ~isempty(S.Ath)
            S.d = size(S.Ath,2);
        elseif ~isempty(S.pB)
            S.d = size(S.pB,2);
        elseif ~isempty(S.pE)
            S.d = size(S.pE,2);
        elseif ~isempty(S.pF)
            S.d = size(S.pf,2);
        else
            S.d = 0;
        end
    end
    if isempty(S.f)
        S.f = zeros(S.n,1);
    end
    if isempty(S.c)
        S.c = 0;
    end
   
    % some other checks
    %sub_error_checks(S);
    
end

if nargin==2 && isfield(S, 'isParametric') && S.isParametric
	% solve the parametric problem for a particular value of the parameter

	% set solver
	switch S.problem_type
		case 'QP',
			S.solver = MPTOPTIONS.qpsolver;
		case 'LP'
			S.solver = MPTOPTIONS.lpsolver;
		case 'LCP'
			S.solver = MPTOPTIONS.lcpsolver;
		otherwise
			error('%s problems not supported.', S.problem_type);
    end

    % check bounds on parameters (issue #108)
    if isfield(S, 'Ath') && ~isempty(S.Ath)
        if any(S.Ath*param_value-S.bth > MPTOPTIONS.abs_tol)
            % the parameter is out of bounds => problem is infeasible
            R.xopt = zeros(S.n,1);
            R.obj = Inf;
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
    end
    
	% substitute the parameter into constraints
	if ~isempty(S.Q)
		% PLCP problem
		S.q = S.q + S.Q*param_value;
		% unset parametric constraints
		S.Q = [];
	else
		% PLP/PQP
		S.b = S.b + S.pB*param_value;
		S.be = S.be + S.pE*param_value;
		S.f = S.f + S.pF*param_value;
		S.c = S.c + param_value'*S.Y*param_value + S.C*param_value;
		% unset parametric constraints
		S.pB = [];
		S.pE = [];
		S.pF = [];
		S.Y = [];
        S.C = [];
	end
	
	% mark the problem as non-parametric
	S.d = 0;
	S.isParametric = 0;
end

if isempty(S.test)
    S.test = 0;
end

% We can not call parametric solvers directly, i.e. without "Opt" and
% "Polyhedron" classes. For this purpose the function "mpt_solvemp" is
% devoted.
if S.d>0   
    error('mpt_solve: Calls to parametric solvers are not supported. Use "mpt_solvemp".');
end


% since MILP/MIQP problems are not in Opt class, we try to
% find out which type of problem we're dealing with
%
% vartype = column array containing the types of the variables.
% 'C' Continuous variable.
% 'I' Integer variable
% 'B' Binary variable
% 'S' Semi-continuous
% 'N' Semi-integer
if isempty(S.problem_type) 
    if ~isempty(S.M)
        S.problem_type = 'LCP';
    else
        if isempty(S.H)
            % if H field is empty
            % LP, MILP, mpLP, mpMILP
            if isempty(S.vartype)
                S.problem_type = 'LP';
                S.vartype = char('C'*ones(S.n,1));
            elseif length(strmatch('C',S.vartype))<length(S.vartype)
                S.problem_type = 'MILP';
            end
        else
            % QP, MIQP, mpQP, mpMIQP
            if isempty(S.vartype)
                S.problem_type = 'QP';
                S.vartype = char('C'*ones(S.n,1));
            elseif length(strmatch('C',S.vartype))<length(S.vartype)
                S.problem_type = 'MIQP';
            end
            
        end
    end
end
% make sure that vartype is a column vector
S.vartype = S.vartype(:);

% initial conditions are zeros if not specified
% can be used by some of the solvers as the starting point
if isempty(S.x0)
    S.x0 = zeros(S.n,1);
end

if isempty(S.A),
    % some solvers don't like when there are no inequality constraints,
    % becuase they use size of the "A" matrix to deduce number of variables.
    % therefore we define A and b as empty matrices of corresponding size.
    S.A = zeros(0, S.n);
    S.b = zeros(0, 0);
end

% check and remove possibly double-sided inequalities
% if ~isempty(S.b)
%     if S.test
%         % do no use global settings
%         [An, bn, Aen, ben] = mpt_ineq2eq(S.A, S.b, 1e-9);
%     else
%         [An, bn, Aen, ben] = mpt_ineq2eq(S.A, S.b);
%     end
%     if ~isempty(ben)
% %        disp('Inequality constraints contain equalities written as double-sided inequalities.');
% %        disp('Extracting double-sided inequalities.');
%         S.A = An;
%         S.m = size(An,1);
%         S.b = bn;
%         S.Ae = [S.Ae; Aen];
%         S.be = [S.be; ben];
%         S.me = size(S.Ae,1);
%     end
% end

if any(S.lb==Inf) || any(S.ub==-Inf)
    % infeasible by construction
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
    return;
end

% replace NaNs with zero rows
S = sub_fix_nan(S);

% % replace any Inf,-Inf terms with S.options.InfBound
% if S.test
%     S.lb(S.lb == -Inf) = -1e9;
%     S.ub(S.ub == Inf) = 1e9;
% else
%     S.lb( S.lb == -Inf | S.lb < -MPTOPTIONS.infbound ) = -MPTOPTIONS.infbound;
%     S.ub( S.ub == Inf | S.ub > MPTOPTIONS.infbound ) = MPTOPTIONS.infbound;
% end

% % remove Inf, -Inf terms from constraints
% minfb = find(S.b == -Inf);
% S.A(minfb,:) = [];
% S.b(minfb) = [];
% pinfb = find(S.b == Inf);
% S.A(pinfb, :) = [];
% S.b(pinfb) = [];
% S.m = S.m-length(minfb)-length(pinfb);

% % remove zero rows from inequality constraints
% if S.test
%     zr = all(abs(S.A)<1e-12, 2);
% else
%     zr = all(abs(S.A)<MPTOPTIONS.zero_tol, 2);
% end
% S.A(zr,:) = [];
% S.b(zr) = [];
% S.m = S.m-nnz(zr);
% 
% % remove zero rows from equality constraints
% if S.test
%     zre = all(abs(S.Ae)<1e-12, 2);
% else
%     zre = all(abs(S.Ae)<MPTOPTIONS.zero_tol, 2);
% end
% if ~isempty(zre)
%     S.Ae(zre,:) = [];
%     S.be(zre) = [];
%     S.me = S.me-nnz(zre);
% end


% make sure that linear part of the objective function is a column vector
S.f = S.f(:);

if any(strcmpi(S.problem_type,{'LP','QP'}))
    if isempty(S.H)
        H = zeros(S.n);
    else
        H = S.H;
    end
    
    % check rank of Ae
    if S.test
        re = rank(full(S.Ae));
    else
        re = rank(full(S.Ae),MPTOPTIONS.rel_tol);
    end
%     % check consistency
%     if S.test
%         rc = rank(full([S.Ae S.be]));
%     else
%         rc = rank(full([S.Ae S.be]),MPTOPTIONS.rel_tol);
%     end

    if re>=S.n
        % overdetermined system, no degrees of freedom
        R.xopt = S.Ae\S.be;
        
        % check constraints
        eqc = ( norm(S.Ae*R.xopt-S.be,Inf) <= MPTOPTIONS.rel_tol );
        if ~isempty(S.b)
            ineqc = ( S.A*R.xopt <= S.b + MPTOPTIONS.rel_tol );
        else
            ineqc = true;
        end
        
        % multipliers        
        % H*x + f + A'*lambda_ineq + Ae'*lambda_eq = 0        
        R.lambda.ineqlin = zeros(S.m,1);
        R.lambda.eqlin = -transpose(S.Ae)\(H*R.xopt + S.f);
        R.lambda.lower = zeros(S.n,1);
        R.lambda.upper = zeros(S.n,1);
        
        if eqc && all(ineqc)
            % feasible solution
            if S.test
                R.exitflag = 1;
            else
                R.exitflag = MPTOPTIONS.OK;
            end
            R.how = 'ok';
        else
            % infeasible
            if S.test
                R.exitflag = 2;
            else
                R.exitflag = MPTOPTIONS.INFEASIBLE;
            end
            R.how = 'infeasible';
        end
        R.obj = 0.5*R.xopt'*H*R.xopt + S.f'*R.xopt + S.c;
        return
    end
%     % for underdetermined system check linearly dependent rows
%     if re<S.me || rc<S.me
%         
%         % if the right hand side is not linearly dependent, infeasible solution
%         if rc>re
%             R.xopt = zeros(S.n,1);
%             R.obj = 0;
%             R.lambda = [];
%             if S.test
%                 R.exitflag = 2;
%             else
%                 R.exitflag = MPTOPTIONS.INFEASIBLE;
%             end
%             R.how = 'infeasible';
%             return
%         end
% 
%         while S.me ~= re
%             % find linearly dependent rows
%             [~,~,p] = lu(sparse([S.Ae S.be]),'vector');
%             rd = p(re+1:end);
%             
%             % remove linearly dependent rows
%             S.Ae(rd,:) = [];
%             S.be(rd) = [];
%             S.me = S.me-length(rd);
%             
%             if S.test
%                 re = rank(full([S.Ae S.be]));
%             else
%                 re = rank(full([S.Ae S.be]),MPTOPTIONS.rel_tol);
%             end
%         end
%         
%     end
    
    % if there are no inequality constraints and only equality constraints,
    % solve the appropriate Lagrange system
    if S.m==0 && isempty(S.lb) && isempty(S.ub) && S.me>0
        if strcmpi(S.problem_type,'QP')
            xopt = [H S.Ae'; S.Ae zeros(S.me)] \ [-S.f; S.be];
            R.xopt = xopt(1:S.n);
            R.lambda.ineqlin = [];
            R.lambda.eqlin = xopt(S.n+1:end);
            R.lambda.lower = [];
            R.lambda.upper = [];
            if S.test
                R.exitflag = 1;
            else
                R.exitflag = MPTOPTIONS.OK;
            end            
            R.how = 'ok';
        else
            R.xopt = S.Ae\S.be;
            R.lambda.ineqlin = [];
            R.lambda.eqlin = -S.Ae'\S.f;
            R.lambda.lower = [];
            R.lambda.upper = [];
            if S.test
                R.exitflag = 3;
            else
                R.exitflag = MPTOPTIONS.UNBOUNDED;
            end            
            R.how = 'unbounded';
        end
        R.obj = 0.5*R.xopt'*H*R.xopt + S.f'*R.xopt + S.c;
        return
    end
end

% if solver is not specified, find the appropriate one and try again
if isempty(S.solver)
    if S.d>0
        S.solver = MPTOPTIONS.parametricsolver;
    else
        S.solver = MPTOPTIONS.([lower(S.problem_type),'solver']);
    end
end

R = sub_callsolver(S, MPTOPTIONS);

% add constant term
if ~isempty(R.obj)
    R.obj = R.obj + S.c;
end

end

%-------------------------------------------------------------------------------
function R=sub_callsolver(S, MPTOPTIONS)

% we prefer CDD and LCP for small LP problems
if isequal(lower(S.solver(1:3)), 'cdd')
	% CDD
    R = mpt_call_cdd(S);
	return

elseif any(strcmpi(S.solver,{'gurobi','gurobi_mex','gurobimex'}))
    % GUROBI
    R = mpt_call_gurobi(S);
    return

elseif isequal(lower(S.solver),'lcp')
    % LCP solver
	if ~isfield(S, 'routine')
		S.routine = [];
	end
	if ~isfield(S, 'H')
		S.H = zeros(S.n);
	end

    R = mpt_call_lcp(S);
    % routine 0 may fail to satisfy constraints, therefore we try to repeat
    % with routine 1 and 2 (will be slower, but more reliable)
    if ~S.test && S.m~=0
        if R.exitflag~=MPTOPTIONS.OK || any(S.A*R.xopt>S.b+MPTOPTIONS.rel_tol)
            if MPTOPTIONS.modules.solvers.lcp.routine==0
                S.routine=1;
                R = mpt_call_lcp(S);
                if R.exitflag~=MPTOPTIONS.OK || any(S.A*R.xopt>S.b+MPTOPTIONS.rel_tol)
                    S.routine=2;
                    R = mpt_call_lcp(S);
                end
%                 if any(S.A*R.xopt>S.b+MPTOPTIONS.rel_tol)
%                     % violation of constraints
%                     R.exitflag = MPTOPTIONS.INFEASIBLE;
%                     R.how = 'infeasible';
%                 end
            elseif MPTOPTIONS.modules.solvers.lcp.routine==1
                S.routine=2;
                R = mpt_call_lcp(S);
%                 if any(S.A*R.xopt>S.b+MPTOPTIONS.rel_tol)
%                     % violation of constraints
%                     R.exitflag = MPTOPTIONS.INFEASIBLE;
%                     R.how = 'infeasible';
%                 end                
            end
        end
	end
	return
end

% make sure S.vartype is defined for other solvers
if ~isfield(S, 'vartype')
	S.vartype = char('C'*ones(S.n,1));
end
if ~isfield(S, 'H')
	% CPLEX needs it
	S.H = zeros(S.n);
end

if strcmpi(S.solver, 'mosek')
    R = mpt_call_mosek(S);

elseif strcmpi(S.solver,'quadprog')
    % quadprog
    R = mpt_call_quadprog(S);
    
elseif any(strcmpi(S.solver,{'glpk','glpkcc','glpkmex'}))
    % GLPK
    R = mpt_call_glpk(S);

elseif any(strcmpi(S.solver,{'nag','nag.e04naf','nag.e04mcf'}))
    % NAG
    R = mpt_call_nag(S);
    
elseif strcmpi(S.solver,'linprog')
    % linprog
    R = mpt_call_linprog(S);
    
elseif any(strcmpi(S.solver,{'cplex','cplexmex','cplexint'}))
    % CPLEX >9
    R = mpt_call_cplex(S);
    
elseif any(strcmpi(S.solver,{'clp','clpmex','mexclp'}))
    % CLP
    R = mpt_call_clp(S);
    
elseif any(strcmpi(S.solver,{'qpOASES','oases','oasesqp'}))
    % qpOASES
    R = mpt_call_qpoases(S);
    
elseif any(strcmpi(S.solver,{'qpc','qpas','qpip','qpc.qpas','qpc.qpip'}))
    % QPC
    R = mpt_call_qpc(S);

elseif strcmpi(S.solver,'qpspline')
    % QPspline solver
    R = mpt_call_qpspline(S);
    
elseif strcmpi(S.solver,'sedumi')
    % SEDUMI
    R = mpt_call_sedumi(S);
    
else
    % for the solver we do not interface directly, call YALMIP
    try
        R = sub_call_yalmip(S);
    catch
        error('mpt_solve: Yalmip error when calling %s solver!',S.solver);
    end
end

% DISABLED
% if the solution is infeasible and rescue option is on, try another
% solver in the list
% if ~S.test
%     s_list = MPTOPTIONS.solvers_list.([S.problem_type,'solvers']);
%     % get position of the next solver in the list
%     pos = find(strcmpi(S.solver,s_list)) + 1;
%     if isempty(pos)
%         return
%     end
%     if R.exitflag ~= 1 && MPTOPTIONS.rescue && pos<=length(s_list)
%         S.solver = s_list{pos};
%         R = mpt_solve(S);
%     end
% end

end

%-------------------------------------------------------------------------------
function R=sub_call_yalmip(S)


% continuous variables
x = sdpvar(S.n,1);
if ~all(S.vartype=='C')
    % binary variables
    x(S.vartype=='B') = binvar(1);
    % integer variables
    x(S.vartype=='I') = intvar(1);
end

% initialize constraint set
if ~isempty(S.A) && ~isempty(S.b)
    Fi = [ S.A*x <= S.b ];
else
    Fi = [];
end
if ~isempty(S.Ae) && ~isempty(S.be)
    Fe = [ S.Ae*x == S.be ];
else
    Fe = [];
end

% objective function
obj = S.f'*x;
if ~isempty(S.H)
    obj = obj + 0.5*x'*S.H*x;
end

%options = MPTOPTIONS.sdpsettings;
%options.solver = S.solver;
options=sdpsettings('Verbose',0,'warning',0,'solver',S.solver);

solution = solvesdp(Fi+Fe, obj, options);
R.xopt = double(x);
R.obj = double(obj);
if ~isempty(S.A) && ~isempty(S.b)
    R.lambda = dual(Fi);
else
    R.lambda = dual(Fe);
end

if solution.problem==0,
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
else
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
end

end

%-------------------------------------------------------------------------------
function S = sub_initialize_args(arg_set, opt_arg_set)
% if the optional argument is not defined, assign empty value

args = [arg_set, opt_arg_set];
S = cell2struct(cell(1,length(args)),args,2);

end

%-------------------------------------------------------------------------------
function sub_error_checks(S)
% additional checks

if isempty(S.M)
    if max(size(S.f))~=S.n
        error('mpt_solve: Field "f" must be a column vector of dimension %d x 1.',S.n);
    end
    validate_realvector(S.f);
    if ~isempty(S.A)
        if size(S.A,2)~=S.n
            error('mpt_solve: Field "A" must be a matrix of dimension %d x %d.',S.m,S.n);
        end
    end
    validate_realmatrix(S.A);
    if ~isempty(S.H)
        if max(size(S.H))~=S.n
            error('mpt_solve: Field "H" must be a square matrix of dimension %d x %d.',S.n);
        end
    end
    validate_realmatrix(S.H);
    if ~isempty(S.Ae)
        if isempty(S.be)
            error('mpt_solve: Field "be" must be provided because "Ae" is given.');
        end
        if size(S.Ae,2)~=S.n
            error('mpt_solve: Field "Ae" must be a matrix of dimension %d x %d.',S.me,S.n);
        end
    end
    validate_realvector(S.c)
    if ~isempty(S.c)
        if ~isscalar(S.c)
            error('mpt_solve: Field "c" must be scalar.');
        end
    end    
    validate_realmatrix(S.Ae);
    if ~isempty(S.be)
        if isempty(S.Ae)
            error('mpt_solve: Field "Ae" must be provided because "be" is given.');
        end
        if size(S.be,2)~=1
            error('mpt_solve: Field "be" must be a column vector of dimension %d x 1,',S.me);
        end
        if size(S.be,1)~=S.me
            error('mpt_solve: Field "be" must be a column vector of dimension %d x 1,',S.me);
        end
    end
    if size(S.b,1)~=S.m
        error('mpt_solve: Field "b" must be a column vector of dimension %d x 1,',S.m);
    end
    validate_realmatrix(S.b);
    if ~isempty(S.lb)
        if size(S.lb,1)~=S.n
            error('mpt_solve: Field "lb" must be a column vector of dimension %d x 1,',S.n);
        end
    end
    validate_realvector(S.lb);
    if ~isempty(S.ub)
        if size(S.ub,1)~=S.n
            error('mpt_solve: Field "ub" must be a column vector of dimension %d x 1,',S.n);
        end
    end
    validate_realvector(S.ub);
    if ~isempty(S.Ath)
        if size(S.Ath,2)~=S.d
            error('mpt_solve: "Ath" must have %d columns.',S.d);        
        end
    end
    validate_realmatrix(S.Ath);
    if size(S.Ath,1)~=size(S.bth,1)
        error('mpt_solve: "Ath" must have same number of rows as "bth".');
    end
    validate_realvector(S.bth);
    if ~isempty(S.pB) && any(size(S.pB)~=[S.m S.d])
        error('mpt_solve: Field "pB" must be of dimension %d x %d,',S.m,S.d);
    end
    validate_realmatrix(S.pB);
    if ~isempty(S.pE) && any(size(S.pE)~=[S.me S.d])
        error('mpt_solve: Field "pE" must be of dimension %d x %d,',S.me,S.d);
    end
    validate_realmatrix(S.pE);
    if ~isempty(S.pF) && any(size(S.pF)~=[S.n S.d])
        error('mpt_solve: Field "pF" must be of dimension %d x %d,',S.n,S.d);
    end
    validate_realmatrix(S.pF);
else
   if size(S.M,1)~=size(S.M,2) || size(S.M,1)~=S.n
       error('mpt_solve: Matrix "M" must be square of dimension %d.',S.n);
   end
   validate_realmatrix(S.M);
   if size(S.q,1)~=S.n
       error('mpt_solve: Field "q" must be a column vector of dimension %d x 1.',S.n);
   end
   validate_realvector(S.q);
end

end

%-------------------------------------------------------------------------------
function S = sub_fix_nan(S)
% replace NaN rows in constraints with zeros

if ~isempty(S.b) || ~isempty(S.A)
    mnan = isnan(S.b) | any(isnan(S.A),2);
    if any(mnan)
        nmnan = nnz(mnan);
        S.A(mnan,:) = zeros(nmnan,S.n);
        S.b(mnan) = zeros(nmnan,1);
    end
end
if ~isempty(S.be) || ~isempty(S.Ae)
    menan = isnan(S.be) | any(isnan(S.Ae),2);
    if any(menan)
        nmenan = nnz(menan);
        S.Ae(menan,:) = zeros(nmenan,S.n);
        S.be(menan) = zeros(nmenan,1);
    end
end

end
