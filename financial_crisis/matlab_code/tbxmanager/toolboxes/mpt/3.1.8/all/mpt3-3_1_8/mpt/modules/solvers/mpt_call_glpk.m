function R = mpt_call_glpk(S)
%
%  MPT_CALL_GLPK: A gateway function to GLPK solver (without errorchecks) 
%  =======================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_glpk(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements call to GLPK solver based on formulation from Opt
%  class. Only LP, MILP and problems are supported. It is assumed that LP/MILP
%  entering this function is of the form 
%                                       T                       
%                                min   f x                  (1) 
%                                                               
%                               s.t.   lb <= x <= ub        (2) 
%                                      Ax <= b              (3) 
%                                                               
%                                      A x = b              (4) 
%                                       e     e                 
%                                      x in {C, I, B }      (5) 
%     where the set {C, I, B}  represents 
%    
%     - C - continuous variables, x in (-oo,oo)   
%     - I - integer variables x in (..., -1, 0, 1, ...)   
%     - B - binary variables x in {0,1}  
%    which is given by strings in vartype field. GLPK accepts this format directly.
%  
%  
%  INPUT
%  -----
%     
%        
%          S              Structure of the Opt class.              
%                         Class: struct                            
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
%          S.vartype      A string specifying the type of          
%                         variable. Supported characters are C     
%                         (continuous), I (integer), B (binary).   
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

% glpk
if ~any(strcmpi(S.problem_type,{'LP','MILP'}))
    error('mpt_call_glpk: GLPK solver does not solve %s problems!',S.problem_type);
end

% options to be passed to glpk (as given in "glpk.m")
%  param = A structure containing the following parameters used to define the
%           behavior of solver.  Missing elements in the structure take on default
%           values, so you only need to set the elements that you wish to change
%           from the default.
%   
%           Integer parameters:
%             msglev (default: 1) 
%                    Level of messages output by solver routines:
%                     0 - No output.
%                     1 - Error messages only.
%                     2 - Normal output.
%                     3 - Full output (includes informational messages).
%   
%             scale (default: 1). Scaling option: 
%                     0 - No scaling.
%                     1 - Equilibration scaling.
%                     2 - Geometric mean scaling, then equilibration scaling.
%                     3 - Geometric then Equilibrium scaling 
%                     4 - Round to nearest power of 2 scaling
%  
%             dual (default: 0). Dual simplex option:
%                     0 - Do not use the dual simplex.
%                     1 - If initial basic solution is dual feasible, use
%                         the dual simplex.
%                     2- Use two phase dual simplex, or if primal simplex 
%                         if dual fails
%  
%             price (default: 1). Pricing option (for both primal and dual simplex):
%                     0 - Textbook pricing.
%                     1 - Steepest edge pricing.
%     
%             r_test (default: 1). Ratio test Technique:
%                     0 - stardard (textbook)
%                     1 - Harris's two-pass ratio test
%     
%             round (default: 0). Solution rounding option:
%   
%                     0 - Report all primal and dual values "as is".
%                     1 - Replace tiny primal and dual values by exact zero.
%  
%             itlim (default: -1). Simplex iterations limit.  
%                   If this value is positive, it is decreased by one each
%                   time when one simplex iteration has been performed, and
%                   reaching zero value signals the solver to stop the search. 
%                   Negative value means no iterations limit.
%   
%             itcnt (default: 200). Output frequency, in iterations.  
%                   This parameter specifies how frequently the solver sends 
%                   information about the solution to the standard output.
%  
%             presol (default: 1). If this flag is set, the routine 
%                    lpx_simplex solves the problem using the built-in LP presolver. 
%                    Otherwise the LP presolver is not used.
%  
%             lpsolver (default: 1) Select which solver to use.
%                    If the problem is a MIP problem this flag will be ignored.
%                     1 - Revised simplex method.
%                     2 - Interior point method.
%                     3 - Simplex method with exact arithmatic.                       
%   
%             branch (default: 2). Branching heuristic option (for MIP only):
%                     0 - Branch on the first variable.
%                     1 - Branch on the last variable.
%                     2 - Branch on the most fractional variable.
%                     3 - Branch using a heuristic by Driebeck and Tomlin.
%  
%             btrack (default: 2). Backtracking heuristic option (for MIP only):
%                     0 - Depth first search.
%                     1 - Breadth first search.
%                     2 - best local bound
%                     3 - Backtrack using the best projection heuristic.
%  
%             pprocess (default: 2) Pre-processing technique option ( for MIP only ):
%                     0 - disable preprocessing
%                     1 - perform preprocessing for root only
%                     2 - perform preprocessing for all levels
%          
%             usecuts (default: 1). ( for MIP only ):
%                    glp_intopt generates and adds cutting planes to
%                    the MIP problem in order to improve its LP relaxation
%                    before applying the branch&bound method 
%                     0 -> all cuts off
%                     1 -> Gomoy's mixed integer cuts
%                     2 -> Mixed integer rounding cuts
%                     3 -> Mixed cover cuts
%                     4 -> Clique cuts
%                     5 -> all cuts
%  
%             binarize (default: 0 ) Binarizeation option ( for mip only ):
%                 ( used only if presolver is enabled )
%                     0 -> do not use binarization
%                     1 -> replace general integer variables by binary ones
%   
%             save (default: 0). If this parameter is nonzero save a copy of 
%                  the original problem to file. You can specify the 
%                  file name and format by using the 'savefilename' and 'savefiletype' 
%                  parameters (see in String Parameters Section here below).
%                  If previous parameters are not defined the original problem 
%                  is saved with CPLEX LP format in the default file "outpb.lp".
%  
%             mpsinfo (default: 1) If this is set, 
%                     the interface writes to file several comment cards, 
%                     which contains some information about the problem. 
%                     Otherwise the routine writes no comment cards. 
%  
%             mpsobj ( default: 2) This parameter tells the 
%                    routine how to output the objective function row: 
%                      0 - never output objective function row 
%                      1 - always output objective function row 
%                      2 - output objective function row if the problem has 
%                          no free rows 
%  
%             mpsorig (default: 0) If this is set, the 
%                     routine uses the original symbolic names of rows and 
%                     columns. Otherwise the routine generates plain names 
%                     using ordinal numbers of rows and columns.
%  
%             mpswide (default: 1) If this is set, the 
%                     routine uses all data fields. Otherwise the routine 
%                     keeps fields 5 and 6 empty. 
%  
%             mpsfree (default: 0) If this is set, the routine 
%                     omits column and vector names every time when possible 
%                     (free style). Otherwise the routine never omits these 
%                     names (pedantic style). 
%  
%   
%           Real parameters:
%             relax (default: 0.07). Relaxation parameter used 
%                   in the ratio test. If it is zero, the textbook ratio test 
%                   is used. If it is non-zero (should be positive), Harris'
%                   two-pass ratio test is used. In the latter case on the 
%                   first pass of the ratio test basic variables (in the case 
%                   of primal simplex) or reduced costs of non-basic variables 
%                   (in the case of dual simplex) are allowed to slightly violate 
%                   their bounds, but not more than relax*tolbnd or relax*toldj 
%                   (thus, relax is a percentage of tolbnd or toldj).
%   
%             tolbnd (default: 10e-7). Relative tolerance used 
%                    to check ifthe current basic solution is primal feasible.
%                    It is not recommended that you change this parameter 
%                    unless you have a detailed understanding of its purpose.
%   
%             toldj (default: 10e-7). Absolute tolerance used to 
%                   check if the current basic solution is dual feasible.  It 
%                   is not recommended that you change this parameter unless 
%                   you have a detailed understanding of its purpose.
%   
%             tolpiv (default: 10e-9). Relative tolerance used 
%                    to choose eligible pivotal elements of the simplex table.
%                    It is not recommended that you change this parameter 
%                    unless you have a detailed understanding of its purpose.
%   
%             objll ( default: -DBL_MAX). Lower limit of the 
%                   objective function. If on the phase II the objective
%                   function reaches this limit and continues decreasing, the
%                   solver stops the search. This parameter is used in the 
%                   dual simplex method only.
%   
%             objul (default: +DBL_MAX). Upper limit of the 
%                   objective function. If on the phase II the objective
%                   function reaches this limit and continues increasing, 
%                   the solver stops the search. This parameter is used in 
%                   the dual simplex only.
%   
%             tmlim (default: -1.0). Searching time limit, in 
%                   seconds. If this value is positive, it is decreased each 
%                   time when one simplex iteration has been performed by the
%                   amount of time spent for the iteration, and reaching zero 
%                   value signals the solver to stop the search. Negative 
%                   value means no time limit.
%   
%             outdly (default: 0.0). Output delay, in seconds. 
%                    This parameter specifies how long the solver should 
%                    delay sending information about the solution to the standard
%                    output. Non-positive value means no delay.
%   
%             tolint (default: 10e-5). Relative tolerance used 
%                    to check if the current basic solution is integer
%                    feasible. It is not recommended that you change this 
%                    parameter unless you have a detailed understanding of 
%                    its purpose.
%   
%             tolobj (default: 10e-7). Relative tolerance used 
%                    to check if the value of the objective function is not 
%                    better than in the best known integer feasible solution.  
%                    It is not recommended that you change this parameter 
%                    unless you have a detailed understanding of its purpose.
%  
%             mipgap (default: 0.0) The relative mip gap tolerance.  If the 
%                    relative mip gap for currently known best integer feasible 
%                    solution falls below this tolerance, the solver terminates 
%                    the search.  This allows obtaining suboptimal interger 
%                    feasible solutions if solving the problem to optimality 
%                    takes too long.
%  
%           String Parameters:
%             savefilename (default: "outpb"). Specify the name to use to 
%                          save the original problem. MEX interface looks for 
%                          this parameter if 'save' parameter is set to 1. If 
%                          no name is provided "outpb" will be used.
%             savefiletype (default: CPLEX format). Specify the format type
%                          used to save the file. Only the following options
%                          are allowed:
%                            'fixedmps' - fixed MPS format (.mps).
%                            'freemps'  - free MPS format (.mps). 
%                            'cplex'    - CPLEX LP format (.lp).
%                            'plain'    - plain text (.txt).
%                            

if ~S.test
    param = MPTOPTIONS.modules.solvers.glpk;
else
    param.msglev=1;
    param.itlim=5000;
    param.lpsolver=1;
    param.scale = 1;
    param.dual = 0;
    param.price = 1;
    param.relax = 0.07;
    param.branch = 2;
    param.btrack = 2;
end

% merge constraints
A = [S.Ae; S.A];
b = [S.be; S.b];
% fixed are first, Ae*x = be, upper-bounded variables are second in A*x<=b,
ctype = char( ['S'*ones(S.me,1); 'U'*ones(S.m,1)] );

% merge lb/ub with inequality constraints
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
    Alb(:,~ilb) = -speye(nnz(~ilb));
    A = [A; Alb];
    b = [b; -S.lb(~ilb)];
    ctype = [ctype; char('U'*ones(nnz(~ilb),1))];
end
if any(~iub)
    Aub = zeros(nnz(~iub),S.n);
    Aub(:,~iub) = speye(nnz(~iub));
    A = [A; Aub];
    b = [b; S.ub(~iub)];
    ctype = [ctype; char('U'*ones(nnz(~iub),1))]; 
end
    
    
% ctype = An array of characters containing the sense of each constraint in the
% constraint matrix.  Each element of the array may be one of the
% following values
% 'F' Free (unbounded) variable (the constraint is ignored).
% 'U' Variable with upper bound ( A(i,:)*x <= b(i)).
% 'S' Fixed Variable (A(i,:)*x = b(i)).
% 'L' Variable with lower bound (A(i,:)*x >= b(i)).
% 'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).

% set LB/UB as Inf because glpkcc accepts Inf values
lb = -Inf(S.n,1);
ub = Inf(S.n,1);
if strcmpi(S.solver,'glpkmex')
    % not tested !!!
    [R.xopt,R.obj,status,lambda_struct] = glpkmex(1,S.f,A,b,ctype,...
         lb,ub,S.vartype,param,method,0);

else 
    % glpkcc
    [R.xopt,R.obj,status,lambda_struct] = glpkcc(S.f,A,b,...
        lb,ub,ctype,S.vartype,1,param);

end

if isempty(lambda_struct.lambda)
    % MILP case - no lambdas
    R.lambda.ineqlin = [];
    R.lambda.eqlin = [];
    R.lambda.lower = [];
    R.lambda.upper = [];
else
    R.lambda.ineqlin = -lambda_struct.lambda(S.me+1:S.me+S.m);
    R.lambda.eqlin = -lambda_struct.lambda(1:S.me);
    if ~isempty(S.lb)
        R.lambda.lower = zeros(S.n,1);
        R.lambda.lower(kept_rows.lb) = -lambda_struct.lambda(S.me+S.m+1:S.me+S.m+numel(kept_rows.lb));
    else
        R.lambda.lower = zeros(S.n,1);
    end
    if ~isempty(S.ub) && isempty(S.lb)
        R.lambda.upper = zeros(S.n,1);
        R.lambda.upper(kept_rows.ub) = -lambda_struct.lambda(S.me+S.m+1:S.me+S.m+numel(kept_rows.ub));
    elseif ~isempty(S.ub) && ~isempty(S.lb)
        R.lambda.upper = zeros(S.n,1);
        R.lambda.upper(kept_rows.ub) = -lambda_struct.lambda(S.me+S.m+numel(kept_rows.lb)+1:S.me+S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
    else
        R.lambda.upper = zeros(S.n,1);
    end
end

switch status
    case {5}   % optimal, feasible
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case {3,4,110,213}   % infeasible
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end
    case {6,214}   % unbounded
        R.how = 'unbounded';
        if S.test
            R.exitflag = 3;
        else
            R.exitflag = MPTOPTIONS.UNBOUNDED;
        end
    case 207
        R.how = ['Maximum of',param.itlim,'reached.'];
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
    case {210,211,212}
        R.how = 'Numerical problems.';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
    case {1,170}
        R.how = 'Other identified error.';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
    otherwise
        R.how = 'infeasible';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
end
