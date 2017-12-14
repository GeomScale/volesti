function R = mpt_call_mosek(S)
% Interface to the MOSEK solver (QP, LP, MILP, MIQP)

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

assert(~strcmpi(S.problem_type, 'LCP'), 'mpt_call_mosek: MOSEK does not solve %s problems!', S.problem_type);

% blc <= a*x <= buc
prob.a = sparse([S.A; S.Ae]);
prob.buc = [S.b; S.be];
prob.blc = [-Inf(S.m, 1); S.be];

% blx <= x <= buc
prob.blx = S.lb;
prob.bux = S.ub;

% QP: min 0.5*x'*q*x + c'*x
switch S.problem_type
    case {'QP', 'MIQP'}
        [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(sparse(S.H)));
        prob.c = S.f(:);
    case {'LP', 'MILP'}
        prob.c = S.f(:);
    otherwise
        error('mpt_cll_mosek: MOSEK does not solve %s problems!', S.problem_type);
end

% integer variables
if (isequal(S.problem_type, 'MILP') || isequal(S.problem_type, 'MIQP')) && ...
        ~isempty(S.vartype)
    % bound binary variables to [0, 1]
    prob.blx(S.vartype~='C') = 0;
    prob.bux(S.vartype=='B') = 1;
    prob.ints.subs = find(S.vartype~='C');
end

% call the solver
command = 'minimize echo(0)'; % keep silent
if S.test
    options = [];
else
    options = MPTOPTIONS.modules.solvers.mosek;
end
[~, res] = mosekopt(command, prob, options);

if res.rcode==2000 || res.rcode==2001
    % infeasible according to callmosek.m
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
    R.how = 'infeasible';
    R.obj = Inf;
    R.xopt = NaN(S.n, 1);
    R.lambda.ineqlin = NaN(S.m, 1);
    R.lambda.eqlin = NaN(S.me, 1);
    R.lambda.lower = NaN(S.n, 1);
    R.lambda.upper = NaN(S.n, 1);
    return
end

if ~isfield(res, 'sol')
    error(res.rmsg);
elseif isfield(res.sol, 'int')
    % integer solution
    out = res.sol.int;
elseif isfield(res.sol, 'bas')
    % prefer the simplex algorithm over interior point
    out = res.sol.bas;
elseif isfield(res.sol, 'itr')
    out = res.sol.itr;
else
    error('mpt_call_mosek: unexpected output from the solver.');
end

switch out.prosta
    case {'PRIMAL_AND_DUAL_FEASIBLE', 'PRIMAL_FEASIBLE'}
        exitflag = 1;
    case 'DUAL_INFEASIBLE'
        exitflag = 3;
    case {'PRIMAL_INFEASIBLE', 'PRIMAL_INFEASIBLE_OR_UNBOUNDED'}
        exitflag = 2;
    case 'MSK_RES_TRM_USER_CALLBACK'
        exitflag = -1;
    case 'MSK_RES_TRM_STALL'
        exitflag = -1;
    case 'UNKNOWN'
        exitflag = -1;
    otherwise
        exitflag = -1;
end
if S.test
    R.exitflag = exitflag;
    R.how = lower(out.solsta);
else
    % translate exitflags to MPT codes
    switch exitflag
        case 1,
            R.exitflag = MPTOPTIONS.OK;
            R.how = 'ok';
        case 2
            R.exitflag = MPTOPTIONS.INFEASIBLE;
            R.how = 'infeasible';
        case 3
            R.exitflag = MPTOPTIONS.UNBOUNDED;
            R.how = 'unbounded';
        otherwise
            R.exitflag = MPTOPTIONS.ERROR;
            R.how = 'error';
    end
end

% optimizer
R.xopt = out.xx(1:S.n);

% objective value
R.obj = out.pobjval;

% Lagrange multipliers
R.lambda.ineqlin = -out.y(1:S.m);
R.lambda.eqlin = -out.y(S.m+1:end);
R.lambda.lower = out.slx;
R.lambda.upper = out.sux;

end
