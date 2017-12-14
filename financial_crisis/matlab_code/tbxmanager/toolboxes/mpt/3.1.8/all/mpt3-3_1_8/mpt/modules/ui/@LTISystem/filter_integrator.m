function filter = filter_integrator(varargin)
% Adds a stedy-state-based integrator for LTI systems

% set up the filter
filter = FilterSetup;
filter.addField('value', true, @islogical);

% the filter impacts the following calls:
filter.callback('addFilter') = @on_addFilter;
filter.callback('constraints') = @on_constraints;

end

%------------------------------------------------
function out = on_addFilter(obj, varargin)
% called when the filter is added

x_ref = obj.x.hasFilter('reference') && isequal(obj.x.reference, 'free');
y_ref = obj.y.hasFilter('reference') && isequal(obj.y.reference, 'free');
u_ref = obj.u.hasFilter('reference') && isequal(obj.u.reference, 'free');

if u_ref
	error('Please remove the "reference" filter from "u" first.');
end
if x_ref && y_ref
	error('Please remove either "x.reference" or "y.reference" filters.');
elseif ~(x_ref || y_ref)
	error('Either state or output tracking must be enabled.');
end

if y_ref
	% we need to enable symbolic references for the state
	if ~obj.x.hasFilter('reference')
		obj.x.with('reference');
	end
	obj.x.reference = 'symbolic';
end

% in either case we need a symbolic reference for the control action
if ~obj.u.hasFilter('reference')
	obj.u.with('reference');
end
obj.u.reference = 'symbolic';

out = [];

end


%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constraints are constructed

out = [];
if ~obj.integrator
	% integrator is not enabled
	return

elseif isequal(obj.x.reference, 'symbolic')
	% output tracking, we need to compute steady-state values for states
	% and inputs from
	%   x_ss = A*x_ss + B*u_ss + f
	%  y_ref = C*x_ss + D*u_ss + g
	
	x_ss = obj.x.Internal.reference.var;
	u_ss = obj.u.Internal.reference.var;
	y_ref = obj.y.Internal.reference.var;
	
	out = out + [ x_ss == obj.A*x_ss + obj.B*u_ss + obj.f ];
	out = out + [ y_ref == obj.C*x_ss + obj.D*u_ss + obj.g ];
	
else
	% state tracking, only compute steady-state inputs from
	%  x_ref = A*x_ref + B*u_ss + f
	x_ref = obj.x.Internal.reference.var;
	u_ss = obj.u.Internal.reference.var;
	
	out = out + [ x_ref == obj.A*x_ref + obj.B*u_ss + obj.f ];
end

end
