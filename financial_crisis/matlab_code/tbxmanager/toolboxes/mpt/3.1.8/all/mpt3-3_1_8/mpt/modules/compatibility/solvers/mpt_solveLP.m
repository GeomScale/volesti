function [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,Aeq,Beq,x0,lpsolver,lb,ub)
%MPT_SOLVELP Interface to various LP solvers
%
% [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,Aeq,Beq,x0,lpsolver,lb,ub)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves an LP problem:
%
%     min  f'x
%     s.t.   A*x <= B
%          Aeq*x = Beq
%          lb <= x <= ub
%
% by using the method specified in 'solver'
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% f        - Optimization objective
% A,B      - Matrices defining inequality constraints
% Aeq,Beq  - Matrices defining equality constraints (optional)
% x0       - Initial value                          (optional)
% lpsolver - Which LP solver to use:
%              lpsolver=0: uses NAG (E04NAG.M)
%              lpsolver=9: uses NAG (E04MBF.M)
%              lpsolver=1: uses linprog.m
%              lpsolver=2: uses CPLEX9 (cplexint)
%              lpsolver=3: uses CDD Criss-Cross (cddmex)
%              lpsolver=4: uses GLPK (glpkmex)
%              lpsolver=5: uses CDD Dual Simplex (cddmex)
%              lpsolver=6: uses SeDuMi
%              lpsolver=7: uses QSopt (qsoptmex)
%              lpsolver=8: uses CPLEX8 (lp_cplex)
%              lpsolver=10: uses XPRESS
%              lpsolver=11: uses MOSEK
%              lpsolver=12: uses OOQP
%              lpsolver=13: uses CLP
%              lpsolver=14: uses BPMPD (bpmpd_mex)
%              lpsolver=15: uses CPLEX (cplexmex)
%              lpsolver=16: uses PDCO
%
% lb, ub   - lower and upper bounds on variables (optional)
%            Note! set "lpsolver=[]" to select the default solver.
%
% Note: if 'lpsolver' is not specified, mptOptions.lpsolver will be used instead
%       (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% xopt      - The optimizer
% fval      - Value of the objective
% lambda    - Vector of Lagrangian multipliers
% exitflag  - An integer specifying result of the optimization:
%                1 - feasible optimal solution found
%                0 - unbounded or undecided problem
%               -1 - infeasible problem
% how       - States the result of optimization ('ok', 'unbounded', 'infeasible')
%
% see also MPT_SOLVEQP, MPT_MPLP, MPT_MPQP

% Copyright is with the following author(s):
%
%(C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch
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

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS=mptopt;
end

if MPTOPTIONS.modules.compatibility.use_old_solveLP
	[xopt,fval,lam,exitflag,how]=mpt_solveLP_26(f,A,B,Aeq,Beq,[],3);
    % equalities are appended at the end
    lambda.ineqlin  = lam(1:size(A,1));
    lambda.eqlin = lam(size(A,1)+1:end);
	return
end

% if true, we bypass validity checks in the Opt class
skip_validation = true;

if skip_validation
	prob = struct('f', f, 'A', A, 'b', B, 'Ae', Aeq, 'be', Beq);
	sol = mpt_solve(prob);
else
	o = Opt('f', f, 'A', A, 'b', B, 'Ae', Aeq, 'be', Beq);
	sol = o.solve();
end
xopt = sol.xopt;

% as of May 21, 2012 the @Opt class does not return correct lagrange
% multipliers because zero rows are removed from the constraints
hack = false;
if hack
	active = find([Aeq; A]*xopt-[Beq; B]>-1e-8);
	lambda = zeros(size([Aeq; A], 1), 1);
	lambda(active)=1;
else
	lambda = sol.lambda;
end

how = sol.how;
exitflag = sol.exitflag;
fval = sol.obj;
if exitflag == MPTOPTIONS.OK
	exitflag = 1;
else
	exitflag = -1;
end
