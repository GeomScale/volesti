function Matrices = mpt_yalmip2mpt(constraints, objective, parameters, requested)
% Returns matrices of the parametric problem
%
%   Matrices = mpt_yalmip2mpt(constraints, objective, parameters, requested)
%
% returns the matrices of the pQP/pLP formulation of the MPC
% optimization problem. If Matrices.qp=1, the problem is
% formulated as
%   min_U 1/2 U'*H*U + x'*F*U + Cf*U + x'*Y*x + Cx*x + Cc
%    s.t. G*U <= W + E*x
%
% If the problem is a pLP (identified by Matrices.qp=0), then:
%   min_U H*U + Cx*x + Cc
%    s.t. G*U <= W + E*x
%
% Here, U is the open-loop sequence of optimal control inputs
% and "x" is the vector of parameters.

% expand the YALMIP model
yopts = sdpsettings('verbose', 0);
w = warning;
warning('off');
try
    [Fexp, failure, cause] = expandmodel(constraints, objective, yopts);
catch
    warning(w);
    error(lasterr);
end
warning(w)
if failure,
    fprintf('\n%s\n\n', cause);
    error('Cannot deal with given setup, see message above.');
end

% now export the expanded model into MPT format
yopts.expand = 0;
[~, ~, ~, model] = export(Fexp, objective, yopts);

% define which variables should be treated as parameters
model.parametric_variables = find(ismember(model.used_variables,getvariables(parameters(:))));
model.requested_variables = find(ismember(model.used_variables,getvariables(requested(:))));

% convert the YALMIP model to MPT matrices
Matrices = yalmip2mpt(model);

% we do not support binaries yet
if ~isempty(model.binary_variables) || ~isempty(model.integer_variables)
    error('Formulations with binary/integer variables not supported.');
end

% Remove equality constraints and trivial stuff from big-M
[equalities,redundant] = mpt_detect_fixed_rows(Matrices);
if ~isempty(equalities)
    % Constraint of the type Ex == W, i.e. lower-dimensional
    % parametric space
    if any(sum(abs(Matrices.G(equalities,:)),2)==0)
        warning('Lower-dimensional constraints in the parametric space.');
        return
    end
end
Matrices = mpt_collect_equalities(Matrices,equalities);
Matrices = mpt_remove_equalities(Matrices,redundant);
Matrices = mpt_project_on_equality(Matrices);

% remove from Matrices fields introduced by YALMIP which might confuse Opt
fields_to_remove = {'requested_variables', 'param_var', 'free_var', ...
    'binary_var_index', 'getback', 'lb', 'ub'};
for i = 1:length(fields_to_remove)
    if isfield(Matrices, fields_to_remove{i})
        Matrices = rmfield(Matrices, fields_to_remove{i});
    end
end
end
