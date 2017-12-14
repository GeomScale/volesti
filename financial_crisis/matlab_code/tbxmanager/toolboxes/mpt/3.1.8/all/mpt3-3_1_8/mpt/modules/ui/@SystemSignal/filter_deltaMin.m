function filter = filter_deltaMin(obj)
%
%  DELTAMIN: Lower bound on the increment of a signal 
%  ===================================================
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
%     This filter introduces the lower bound lb  on the difference of a signal lb
%  <= Delta s_k . The signal difference is defined as Delta s_k = s_k - s_k-1 
%  where s  is the system signal and k  is the time instant. The filter is
%  activated by calling s.with('deltaMin').
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
 
 
filter = FilterSetup;
filter.addField('value', -Inf(obj.n, 1), @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('addFilter') = @on_addFilter;
filter.callback('removeFilter') = @on_removeFilter;
filter.callback('set') = @on_set;

end

%-----------------------------------------------
function obj = on_set(obj, value)
% called when the filter's values are changed

error(validate_vector(value, obj.n, 'value of deltaMin'));
if obj.hasFilter('deltaMax') && any(value > obj.deltaMax)
	error('Lower bound cannot exceed upper bound.');
end
obj.deltaMin = value;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when creating constraints

out = [];
if ~obj.isKind('x')
	% non-state signals have the previous value in
	% obj.Internal.previous.var
	previous = obj.Internal.previous.var;

	% constraints
	delta = obj.var(:, 1)-previous;
	for i = 1:obj.n
		% Do not include +/-Inf bounds
		if ~isinf(obj.deltaMin(i))
			out = out + [ delta(i) >= obj.deltaMin(i) ];
		end
	end
end

for k = 2:obj.N
	delta = obj.var(:, k)-obj.var(:, k-1);
	for i = 1:obj.n
		% Do not include +/-Inf bounds
		if ~isinf(obj.deltaMin(i))
			out = out + [ delta(i) >= obj.deltaMin(i) ];
		end
	end
end	

end

%------------------------------------------------
function out = on_addFilter(obj, varargin)
% called after the filter was added

% non-state signals require the previous value to be stored in vector of
% initial conditions
if ~obj.isKind('x')
	if ~obj.hasFilter('previous')
		obj.addFilter('previous');
	end
	% register ourselves with the "previous" filter (only register
	% once)
	obj.previous = union(obj.previous, mfilename);
end
out = [];

end

%------------------------------------------------
function out = on_removeFilter(obj, varargin)
% called prior to the filter is removed

% non-state signals require the previous value to be stored in vector of
% initial conditions
if ~obj.isKind('x')
	% deregister ourselves from the "previous" filter
	obj.previous = setdiff(obj.previous, mfilename);
	if isempty(obj.previous)
		% we were the last ones requiring the previous value, remove the
		% "init" filter
		obj.removeFilter('previous');
	end
end
out = [];

end
