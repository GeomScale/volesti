function filter = filter_PWApenalty(varargin)
% PWA penalty
%
% Adds a PWA penalty on the k-th step of the prediction horizon. The
% penalty function must be scalar-valued.
%
% This example adds a PWA penalty on the final predicted state, assumes the
% PWA function is given as the function 'obj' of polyunion "PU".
%
%   model.x.with('PWApenalty')
%   model.x.PWApenalty.step = N+1; % N+1 is the index of the terminal state
%   model.x.PWApenalty.polyunion = PU;
%   model.x.PWApenalty.function = 'obj';
%   model.x.PWApenalty.isconvex = true;
%
% Note: this filter conflicts with standard "penalty" filters (either
% 'penalty' or 'terminalPenalty'). It is the user's responsibility to
% switch off the standard filters manually (or set the corresponding
% standard penalties to zeros).

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt;
end

% set up the filter
filter = FilterSetup;
filter.addField('step', []);
filter.addField('polyunion', []);
filter.addField('function', '');
filter.addField('isconvex', []);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('objective') = @on_objective;
filter.callback('instantiate') = @on_instantiate;
filter.callback('uninstantiate') = @on_uninstantiate;
filter.callback('getVariables') = @on_variables;

end

%------------------------------------------------
function out = on_variables(obj, varargin)
% called when filter's variables are requested

% Response: structure (or an array of structures) with following fields:
%
%  .var: sdpvar representation of the introduced variable
%  .parametric: logical, if true, the variable will become part of the
%              vector of initial conditions

% return the variable which represents epigraphs of convex functions
out.var = obj.Internal.PWApenalty_epigraph;
out.parametric = false;

end

%------------------------------------------------
function out = on_instantiate(obj, varargin)
% called after the object was instantiated

% validate arguments
local_validate(obj);

% PWA penalties are modeled using an additional epigraph variable
obj.Internal.PWApenalty_epigraph = sdpvar(1, 1);
out = [];

end

%------------------------------------------------
function out = on_uninstantiate(obj, varargin)
% called when the YALMIP representation of variables is removed

% clear the internal variable
obj.Internal.PWApenalty_epigraph = [];
out = [];

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when creating constraints

if ~obj.PWApenalty.isconvex
	error('Only convex PWA penalties are supported.');
end

e = obj.Internal.PWApenalty_epigraph;
k = obj.PWApenalty.step;
out = [];
% create the epigraph
for i = 1:obj.PWApenalty.polyunion.Num
	% obtain local expressions of the PWA function
	fun = obj.PWApenalty.polyunion.Set(i).getFunction(obj.PWApenalty.function);
	% add onstraints "e >= F_i*x + g_i"
	out = out + [ e >= fun.F*obj.var(:, k) + fun.g ];
end

% require the variable to belong to the domain of the funciton
out = out + [ ismember(obj.var(:, k), obj.PWApenalty.polyunion.convexHull) ];

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when creating the objective function

% just minimize the epigraph variable
out = obj.Internal.PWApenalty_epigraph;
        
end

%------------------------------------------------
function local_validate(obj)
% validate filter's inputs

% check "isconvex"
if ~islogical(obj.PWApenalty.isconvex)
	error('The "isconvex" parameter must be a logical.');
end

% "step" must be an integer
if ~isnumeric(obj.PWApenalty.step) || mod(obj.PWApenalty.step, 1)~=0
	error('The "step" parameter must be an integer.');
end	

% function must be a string
if ~ischar(obj.PWApenalty.function)
	error('The "function" parameter must be a string.');
end

% "polyunion" must be correct
if ~isa(obj.PWApenalty.polyunion, 'PolyUnion')
	error('The "polyunion" parameter must be a PolyUnion object.');
end
if ~obj.PWApenalty.polyunion.hasFunction(obj.PWApenalty.function)
	error('No such function "%s" in the polyunion', obj.PWApenalty.function);
end
fun = obj.PWApenalty.polyunion.Set(1).getFunction(obj.PWApenalty.function);
if ~isa(fun, 'AffFunction')
	error('The function "%s" must be PWA.', obj.PWApenalty.function);
end
if fun.D~=obj.n
	error('Function "%s" must map from R^%d.', obj.PWApenalty.function, obj.n);
end
if fun.R~=1
	error('Function "%s" must be scalar-valued.', obj.PWApenalty.function);
end

end
