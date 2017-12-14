function [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mpqp(Matrices,Options)
%MPT_MPQP Explicitly solves the given quadratic program (QP)
%
% [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mpqp(Matrices,Options)
%
% NOTE: This function is obsolete. Use the @Opt class to define and solve
% mpLP problems. 
%
% This method simply converts the parametric problem into an instance of
% the Opt class, which will solve it using the default parametric solver
% (either PLCP or MPQP).
%
% See "help mpt_mpqp_26" for a detailed description of input/output
% arguments.

% Copyright is with the following author(s):
%
% (C) 2012 Michal Kvasnica, Slovak University of Technology in Bratislava
%     michal.kvasnica@stuba.sk

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

narginchk(1, 2);

% This function acts as a wrapper for the Opt class and delegates solving
% the parametric program to mpt_solvemp().

% make sure all parameters are bounded
if isfield(Matrices,'lb')
    Matrices.lb(Matrices.lb==-Inf) = -MPTOPTIONS.infbound;
end
if isfield(Matrices,'ub')
    Matrices.ub(Matrices.ub==Inf) = MPTOPTIONS.infbound;
end

% remove from Matrices fields introduced by YALMIP which might confuse Opt
fields_to_remove = {'requested_variables', 'param_var', 'free_var', ...
	'binary_var_index', 'getback'};
for i = 1:length(fields_to_remove)
	if isfield(Matrices, fields_to_remove{i})
		Matrices = rmfield(Matrices, fields_to_remove{i});
	end
end

% % exit immediately if problem is infeasible
% %
% % TODO: this check should be part of Opt/validate()
% P = Polyhedron('A', [Matrices.G -Matrices.E], 'b', Matrices.W);
% if ~P.isFullDim
% 	disp('Infeasible optimization problem from the begining');
% 	Pn = polytope;
% 	Fi = {};
% 	Gi = {};
% 	activeConstraints = {};
% 	Phard = polytope;
% 	details = [];
% 	return
% end


% use Opt/solve() to solve the problem
mplp = Opt(Matrices);
solution = mplp.solve();

% Convert the solution to MPT2 format:
% polyhedral partition
if solution.xopt.Num>0
    Pn = polytope(solution.xopt.Set);
else
    % infeasible
    Pn = polytope;
end
    

% active constraints (not returned for now)
activeConstraints = {};

% set of feasible parameters
try
    Phard = polytope(solution.xopt.convexHull);
catch
    Phard = polytope(Opt(Matrices).feasibleSet(solution.xopt.Set));
end

% optimizer
Fi = cell(1,solution.xopt.Num); 
Gi = cell(1,solution.xopt.Num);
if solution.xopt.Num>0
    optimizer = solution.xopt.Set.getFunction('primal');
    for i = 1:length(optimizer)
        Fi{i} = optimizer(i).F;
        Gi{i} = optimizer(i).g;
    end
end

% cost
details.Ai = cell(1,solution.xopt.Num);
details.Bi = cell(1,solution.xopt.Num);
details.Ci = cell(1,solution.xopt.Num);
if solution.xopt.Num>0
    cost = solution.xopt.Set.getFunction('obj');
    for i = 1:length(cost)
        if isa(cost(i), 'QuadFunction')
            details.Ai{i} = cost(i).H;
        end
        details.Bi{i} = cost(i).F;
        details.Ci{i} = cost(i).g;
    end
end

end
