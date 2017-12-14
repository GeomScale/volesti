function filter = filter_previous(varargin)
% Value of a signal in the previous time step

% set up the filter
filter = FilterSetup;
filter.hidden = true;
filter.transient = true; % transient filters are not saved
filter.addField('value', {});

% the filter impacts the following calls:
filter.callback('instantiate') = @on_instantiate;
filter.callback('uninstantiate') = @on_uninstantiate;
filter.callback('constraints') = @on_constraints;
filter.callback('getVariables') = @on_variables;

end

%------------------------------------------------
function out = on_variables(obj, varargin)
% called when filter's variables are requested
%
% Response: structure (or an array of structures) with following fields:
%
%  .var: sdpvar representation of the introduced variable
%  .parametric: logical, if true, the variable will become part of the
%               vector of initial conditions

if isa(obj.Internal.previous.var, 'sdpvar')
	out.var = obj.Internal.previous.var;
	out.parametric = true;
else
	out = [];
end

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

out = [];
if isa(obj.Internal.previous.var, 'sdpvar')
	% bound symbolic variables using signal's min/max values
	%
	% do not include +/-Inf bounds
	v = obj.Internal.previous.var;
	for i = 1:length(v)
		if ~isinf(obj.min(i))
			out = out + [ obj.min(i) <= v(i) ];
		end
		if ~isinf(obj.max(i))
			out = out + [ v(i) <= obj.max(i) ];
		end
	end
end

end

%------------------------------------------------
function out = on_instantiate(obj, varargin)
% called after the object was instantiated

if obj.isKind('x')
	% state signals already have the initial value in obj.var(:, 1)
	obj.Internal.previous.var = [];
else
	% create new variable for other signals
	obj.Internal.previous.var = sdpvar(obj.n, 1);
end
out = [];

end

%------------------------------------------------
function out = on_uninstantiate(obj, varargin)
% called when the YALMIP representation of variables is removed

if isa(obj.Internal.previous.var, 'sdpvar')
	obj.Internal.previous.var = [];
end
out = [];

end
