function filter = filter_softMin(obj)
%
%  SOFTMIN: Soft lower bound constraint 
%  =====================================
%  
%  
%  SYNTAX
%  ------
%     
%      
%    
%  
%  DESCRIPTION
%  -----------
%     This filter will soften the lower bound constraint on a selected signal in
%  the MPC setup. Without this filter, lower bounds are hard, i.e., the signal has
%  to satisfy z_min <= z_k. With this filter added, the lower bound is soft and can
%  be violated by some positive margin d_k, i.e., the new constraint becomes z_min
%  - d_k <= z_k. The slack variables d_k  are then penalized in the MPC cost
%  function by adding the term || Q d_k||_p.
%    To enable this filter, first call model.signal.with('softMin'), where you
%  replace signal by the actual system's signal (typically by x for state
%  variables, u for inputs, and y for outputs).
%    With the filter enabled, you can set the maximal allowed violation of the
%  constraint in model.signal.softMin.maximalViolation, and specify penalization of
%  the slack variables by setting model.signal.softMin.penalty.
%    To remove this filter, call model.signal.without('softMin').
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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
if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt;
end

% set up the filter
filter = FilterSetup;
filter.addField('penalty', AffFunction(MPTOPTIONS.infbound*ones(1, obj.n)), @(x) isa(x, 'Function'));
filter.addField('maximalViolation', 1e3*ones(obj.n, 1), @isnumeric);

% this filter depends on the "min" filter
filter.dependsOn('min') = true;

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('objective') = @on_objective;
filter.callback('instantiate') = @on_instantiate;
filter.callback('uninstantiate') = @on_uninstantiate;
filter.callback('addFilter') = @on_addFilter;
filter.callback('removeFilter') = @on_removeFilter;
filter.callback('getVariables') = @on_variables;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_variables(obj, varargin)
% called when filter's variables are requested

% Response: structure (or an array of structures) with following fields:
%
%  .var: sdpvar representation of the introduced variable
%  .parametric: logical, if true, the variable will become part of the
%              vector of initial conditions
if isa(obj.Internal.soft_min, 'sdpvar')
	out.var = obj.Internal.soft_min;
	out.parametric = false;
else
	out = [];
end

end

%------------------------------------------------
function out = on_instantiate(obj, varargin)
% called after the object was instantiated

% soft constraint require introducing new variables
obj.Internal.soft_min = sdpvar(obj.n, obj.N, 'full');
out = [];

end

%------------------------------------------------
function out = on_uninstantiate(obj, varargin)
% called when the YALMIP representation of variables is removed

% soft constraint require introducing new variables
obj.Internal.soft_min = [];
out = [];

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when creating constraints

s = obj.Internal.soft_min;
out = [];
for i = 1:obj.n
	out = out + [ obj.min(i, :) - s(i, :) <= obj.var(i, :) ];
	out = out + [ 0 <= s(i, :) <= obj.softMin.maximalViolation(i, :) ];
end

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when creating the objective function

s = obj.Internal.soft_min;
out = 0;
for k = 1:obj.N
	out = out + obj.softMin.penalty.feval(s(:, k));
end
        
end


%------------------------------------------------
function obj = on_addFilter(obj)
% called after the filter was added

% we need to deactivate the "min" filter
obj.disableFilter('min');

end

%------------------------------------------------
function obj = on_removeFilter(obj)
% called prior to the filter is removed

% we need to re-activate the "min" filter
obj.enableFilter('min');

end

%------------------------------------------------
function obj = on_set(obj, value)
% validation

error(validate_vector(value.maximalViolation, obj.n, 'maximal violation'));
error(obj.validatePenalty(value.penalty));

obj.softMin.maximalViolation = value.maximalViolation;
obj.softMin.penalty = value.penalty;

end
