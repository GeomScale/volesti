function [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,Aeq,Beq,x0,solver,options,rescue)
%MPT_SOLVEQP Interface to various QP solvers
%
% [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,Aeq,Beq,x0,solver,options);
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves a QP problem:
%
%     min  0.5*x'*H*x+f'x
%     s.t. A x <= B
%          Aeq x = Beq
%
% by using the method specified in SOLVER
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% H,f      - Optimization objective
% A,B      - Matrices defining inequality constraints
% Aeq,Beq  - Matrices defining equality constraints
% x0       - Initial value
% solver   - Which QP solver to use:
%              solver=0:  uses NAG (E04NAF.M)
%              solver=1:  uses QUADPROG.M
%              solver=2:  uses CPLEX 9 (cplexint) 
%              solver=3:  uses SeDuMi
%              solver=4:  uses CPLEX 8 (qp_cplex)
%              solver=5:  uses XPRESS
%              solver=6:  uses MOSEK
%              solver=7:  uses OOQP
%              solver=8:  uses CLP (mexclp)
%              solver=9:  uses BPMPD
%              solver=10: uses CPLEX (cplexmex)
%
% options  - options set by 'optimset' function (only for quadprog)
%
% Note: if 'solver' is not specified, mptOptions.qpsolver will be used instead
%       (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% xopt      - The optimizer
% lambda    - Vector of Lagrangian multipliers
% how       - States the result of optimization ('ok', 'unbounded', 'infeasible')
% objqp     - Value of the objective function at the optimizer
%
% see also MPT_SOLVELP, MPT_MPQP

% Copyright is with the following author(s):
%
%(C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%              kvasnica@control.ee.ethz.ch
%(C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%         baotic@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

narginchk(4, 10);

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<6,
    Aeq = [];
    Beq = [];
end

if nargin<7,
    x0=[];
end
if nargin<8
    solver=mptOptions.qpsolver;
end
if nargin<10,
    rescue = mptOptions.rescueQP;
end

if isempty(A),
    % some solvers don't like when there are no inequality constraints,
    % becuase they use size of the "A" matrix to deduce number of variables.
    % therefore we define A and b as empty matrices of corresponding size.
    nx = length(f);    
    A = zeros(0, nx);
    B = zeros(0, 0);
end

% some solvers don't like +/- Inf terms in constraints
minfb = find(B == -Inf);
if ~isempty(minfb),
    % trivially infeasible problem
    xopt = repmat(NaN, length(f), 1);
    objqp = NaN;
    lambda = [];
    exitflag = -1;
    how = 'infeasible';
    return
end
pinfb = find(B == Inf);
A(pinfb, :) = [];
B(pinfb) = [];

[m,n]=size(A);
[meq,neq]=size(Aeq);
f = f(:);

if rescue,
    H_orig = H;
    f_orig = f;
    A_orig = A;
    B_orig = B;
    Aeq_orig = Aeq;
    Beq_orig = Beq;
    x0_orig = x0;
end

if solver==-1
    % no QP solver available
    fprintf('\n\n')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('There is no QP solver installed on your system.');
    disp('You will not be able to solve problems with quadratic cost function.');
    disp('Please be sure that you set ''probStruct.norm = 1'' before calling any control routine!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n')
    error('mpt_solveQP: No QP solver available');
    
elseif solver==0
    % NAG QP solver
    
    
    % convert sparse matrices to full matrices
    if issparse(H),
        H = full(H);
    end
    if issparse(f)
        f = full(f);
    end
    if issparse(A)
        A = full(A);
    end
    if issparse(B)
        B = full(B);
    end
    if issparse(Aeq)
        Aeq = full(Aeq);
    end
    if issparse(Beq)
        Beq = full(Beq);
    end
    
    [m,n] = size(A);
    itmax=450;
    if isempty(x0),
        x0=zeros(n,1);
    end
    
    bl=[-1e8*ones(m+n,1); Beq];
    bu=[1e8*ones(n,1);B; Beq];
    lp=0;
    cold=1;
    istate = zeros(length(bu),1);
    featol = sqrt(eps)*ones(length(bu),1);
    msglvl=-1;
    bigbnd = 1e10;
    orthog = 1;
    ifail=1;
    [xopt,iter,objqp,clambda,istate,ifail] = e04naf(bl,bu,'mpt_qphess',x0,f(:),[A; Aeq],H,...
        lp,cold,istate,featol,msglvl,itmax,bigbnd,orthog,ifail);
    
    switch ifail
    case 0
        how = 'ok';
        exitflag=1;
    case 2
        how = 'unbounded';
        exitflag=-1;
    case 6
        how = 'infeasible';
        exitflag=-1;
    case {1,3}
        % local minimum
        how = 'ok';
        exitflag=-1;
    case {9}
        error('mpt_solveQP: An input parameter is invalid.');
    otherwise
        how = 'other';
        exitflag=0;
    end
    lambda=-clambda(n+1:m+n);
    
    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

    
elseif solver==1,
    % quadprog
    
    % convert sparse matrices to full matrices
    if issparse(H),
        H = full(H);
    end
    if issparse(f)
        f = full(f);
    end
    if issparse(A)
        A = full(A);
    end
    if issparse(B)
        B = full(B);
    end
    if issparse(Aeq)
        Aeq = full(Aeq);
    end
    if issparse(Beq)
        Beq = full(Beq);
    end
    
    if nargin<9,
        options=optimset(optimset('quadprog'),'Display','off','LargeScale','off', 'Algorithm', 'active-set');
    elseif isempty(options),
        options=optimset(optimset('quadprog'),'Display','off','LargeScale','off', 'Algorithm', 'active-set');
    end
    options.Display = 'off';
    
    % make sure the Hessian is symmetric.
    if norm(H-H', Inf) < 1e-10,
        % but we only remove numerical noise if the hessian is only
        % "slightly" wrong. for clerly non-symmetrical hessians we still
        % let quadprog to display a proper warning
        H = (H + H')*0.5;
    end
    
    [xopt,objqp,exitflag,OUTPUT,lambdav]=quadprog(H,f,A,B,Aeq,Beq,[],[],x0,options);
    if exitflag>0 %then QUADPROG converged with a solution X.
        how = 'ok';
    else
        % ==0 then the maximum number of iterations was exceeded (only occurs
        %     with large-scale method).
        % < 0 then the problem is unbounded, infeasible, or 
        %     QUADPROG failed to converge with a solution X.  else
        how = 'infeasible';
    end
    lambda=lambdav.ineqlin;

    if exitflag <= 0 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

    
elseif solver==2,
    % CPLEX 9
    nc = size(A,1);
    nx = size(A,2);
    A = [A; Aeq];
    B = [B; Beq];
    INDEQ = (nc+1:size(A,1))';  % indices of equality constraints
    LB = [];
    UB = [];
    OPTIONS.verbose = 0;
    OPTIONS.lic_rel = 1e2;  % after how many runs to release the license    
    PARAM.int=[1063, 2];    % use the dual-simplex method
    VARTYPE = []; % all variables are continuous
    
    [xopt,objqp,SOLSTAT,DETAILS] = cplexint(H, f(:), A, B, INDEQ, [], LB, UB, VARTYPE, PARAM, OPTIONS);
    
    exitflag = -1;
    lambda = -DETAILS.dual;
    how = lower(DETAILS.statstring);
    if strcmp(how, 'optimal') | strcmp(how, 'optimalrelaxed') | strcmp(how, 'optimaltol'),
        how = 'ok';
        exitflag = 1;
    end

    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

    
elseif solver==4,
    % CPLEX 8
    ctype= [char('L'*ones(m,1)); char('E'*ones(meq,1))];
    A = [A; Aeq];
    B = [B; Beq];
    OPTIONS.int=[1063, 2];    % use the dual-simplex method
    [xopt,objqp,status,slack,lambda]=qp_cplex(1,H,f(:),A,B,ctype,[],[],[],[],0,OPTIONS);
    
    % Matlab-to-CPLEX interface ver 1.0 used this call:
    % [xopt,objqp,status,lambda] = qp_cplex(1,f(:),H,A,B,ctype,[],[],[],0,0);
    %
    % Note: get the latest version of the interface from:
    % http://control.ee.ethz.ch/~hybrid/cplexint.msql

    switch status
    case {1}   % feasible
        how = 'ok';
        exitflag=1;
    case {0,2} % infeasible
        how = 'infeasible';
        exitflag=-1;
    case {3}   % unbounded
        how = 'unbounded';
        exitflag=-1;
    otherwise
        how = 'infeasible';
        exitflag=-1;
    end

    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

elseif solver==3,
    % SeDuMi / YALMIP
    [xopt,lambda,how,exitflag,objqp]=yalmipQP(H,f,A,B,Aeq,Beq,x0,'sedumi');
    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

elseif solver==5,
    % XPress / YALMIP
    [xopt,lambda,how,exitflag,objqp]=yalmipQP(H,f,A,B,Aeq,Beq,x0,'xpress');
    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

elseif solver==6,
    % Mosek / YALMIP
    [xopt,lambda,how,exitflag,objqp]=yalmipQP(H,f,A,B,Aeq,Beq,x0,'mosek');
    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

elseif solver==7,
    % OOQP / YALMIP
    [xopt,lambda,how,exitflag,objqp]=yalmipQP(H,f,A,B,Aeq,Beq,x0,'ooqp');
    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

elseif solver==8
    % CLP
    [xopt, lambda, status] = clp(H, f, A, B, Aeq, Beq);

    if status==0,
        how = 'ok';
        exitflag = 1;
    elseif status==1,
        how = 'infeasible';
        exitflag = -1;
    else
        how = 'unbounded';
        exitflag = -1;
    end

    objqp = 0.5*xopt'*H*xopt + f'*xopt;

    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

    
elseif solver==9
    % BPMPD
    
    nA = size(A, 1);
    nAeq = size(Aeq, 1);
    An = sparse([A; Aeq]);
    Bn = [B; Beq];
    ind = [repmat(-1, nA, 1); repmat(0, nAeq, 1)];
    
    [xopt,y,s,w,howout] = bp(H, An, full(Bn), full(f(:)), ind, [], [], [], []);
    lambda = -y;

    objqp = 0.5*xopt'*H*xopt + f'*xopt;
    
    switch howout,
        case 'optimal solution'
            exitflag = 1;
        otherwise
            exitflag = -1;
    end
    
    if exitflag==1,
        how = 'ok';
        return
    else
        how = 'infeasible';
    end

    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end

    
elseif solver==10
    % CPLEX interfaced with CPLEXMEX (by Nicolo Giorgetti)
    
    SENSE = 1; % minimize
    F = f(:);  % linear objective must be a column vector

    [nc, nx] = size(A);
    A = [A; Aeq];
    B = [B; Beq];
    
    % type of constraints
    % 'L' - less than or equal (<=)
    % 'E' - equality constraint (=)
    CTYPE = [repmat('L', nc, 1); repmat('E', size(Aeq, 1), 1)];
    
    % all variables are continuous
    VARTYPE = repmat('C', nx, 1);
    
    % lower an upper bounds on optimization variables
    LB = [];
    UB = [];

    PARAM.errmsg = 0;
    PARAM.niter = 1e3; % release CPLEX license after 1000 runs of the QP algorithm

    [xopt,objqp,exitflag,DETAILS] = cplexmex(SENSE,H,F,A,B,CTYPE,LB,UB,VARTYPE,x0,PARAM);

    lambda = -DETAILS.lambda;
    
    if exitflag==1,
        how = 'ok';
        return
    else
        exitflag = -1;
        how = 'infeasible';
    end

    if exitflag ~= 1 & rescue,
        % solution is not optimal, try another solver
        qpsolvers = mptOptions.solvers.qp;
        
        % get position of curent solver
        curpos = find(qpsolvers==solver);
        
        if curpos < length(qpsolvers)
            % if this is not the last solver in the list, try other solver
            nextsolver = qpsolvers(curpos+1);
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H_orig, f_orig, A_orig, B_orig,...
                Aeq_orig, Beq_orig, x0_orig, nextsolver, [], 1);
        end
    end
    
else
    error('mpt_solveQP: unknown value for parameter "solver"!');
end




%-------------------------------------------------------------------------------
function [xopt,lambda,how,exitflag,objqp]=yalmipQP(H,f,A,B,Aeq,Beq,x0,solver)

global mptOptions

[nc,nx] = size(A);
f = f(:);
x = sdpvar(nx,1);
F = set(A*x <= B);
if ~isempty(Aeq),
    F = F + set(Aeq*x == Beq);
end

options = mptOptions.sdpsettings;
options.solver = solver;
%options=sdpsettings('Verbose',0,'warning',0,'solver',solver);

solution = solvesdp(F, 0.5*x'*H*x + f'*x, options);
xopt = double(x);
objqp = 0.5*xopt'*H*xopt + f'*xopt;
lambda = [A;Aeq]*xopt - [B; Beq];

if solution.problem==0,
    how = 'ok';
    exitflag = 1;
else
    how = 'infeasible';
    exitflag = -1;
end
