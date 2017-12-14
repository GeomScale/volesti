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

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS=mptopt;
end

if MPTOPTIONS.modules.compatibility.use_old_solveQP
	[xopt,lam,how,exitflag,objqp]=mpt_solveQP_26(H,f,A,B,Aeq,Beq,x0,solver);
    % equalities are appended at the end
    lambda.ineqlin  = lam(1:size(A,1));
    lambda.eqlin = lam(size(A,1)+1:end);
	return
end

% if true, we bypass validity checks in the Opt class
skip_validation = true;

if skip_validation
	prob = struct('H', H, 'f', f, 'A', A, 'b', B, ...
		'Ae', Aeq, 'be', Beq, 'lb', [], 'ub', [], 'quickqp', true);
	sol = mpt_solve(prob);
else
	o = Opt('H', H, 'f', f, 'A', A, 'b', B, 'Ae', Aeq, 'be', Beq);
	sol = o.solve();
end

if false
	% for debugging
	if isequal(how, 'ok') && norm(xopt-sol.xopt)>1e-8
		how;
	end
	if isequal(how, 'ok') && length(lambda)~=length(sol.lambda)
		how;
	end
	if ~isequal(how, sol.how)
		how;
	end
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
objqp = sol.obj;
if exitflag == MPTOPTIONS.OK
	exitflag = 1;
else
	exitflag = -1;
end
