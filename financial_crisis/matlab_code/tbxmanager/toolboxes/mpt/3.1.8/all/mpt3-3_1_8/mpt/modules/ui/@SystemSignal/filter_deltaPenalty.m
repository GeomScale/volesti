function filter = filter_deltaPenalty(varargin)
%
%  DELTAPENALTY: Penalizes the increment of a signal 
%  ==================================================
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
%  s_k-1 where s  is the system signal and k  is the time instant. The filter is
%  activated by calling s.with('deltaPenalty').
%    Any signal can have this filter, e.g. we can without problems enable slew-rate
%  constraints on states, outputs, or even on binary variables. Non-state signals,
%  however, require that the previous value is specified when calling
%  MPCController/evaluate, e.g. ctrl.evaluate(x0, 'u.previous').
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
filter.addField('value', []);

% the filter impacts the following calls:
filter.callback('objective') = @on_objective;
filter.callback('set') = @on_set;
filter.callback('addFilter') = @on_addFilter;
filter.callback('removeFilter') = @on_removeFilter;

end

%------------------------------------------------
function out = on_addFilter(obj, varargin)

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

%------------------------------------------------
function out = on_objective(obj, varargin)

out = 0;
if isempty(obj.deltaPenalty)
	return
end

if ~obj.isKind('x')
	% non-state signals require a symbolic previous value
	previous = obj.Internal.previous.var;
	% penalization of the first step
	out = out + obj.deltaPenalty.feval(obj.var(:, 1)-previous);
	% now penalize remaining steps
end

% penalize increments
for k = 2:obj.N
	out = out + obj.deltaPenalty.feval(obj.var(:, k) - obj.var(:, k-1));
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called prior to property being set

% validate the penalty (empty penalty means no penalization)
if ~isempty(P)
	error(obj.validatePenalty(P));
end
obj.deltaPenalty = P;

end
