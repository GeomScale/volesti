function filter = filter_terminalPenalty(obj)
%
%  TERMINALPENALTY: Penalizes the final predicted state in the MPC problem 
%  ========================================================================
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
%     This filter adds the term  || Q x_N ||_p  to the MPC cost function.
%  Properties of the penalty (i.e., the weighting matrix Q  and the norm p) are
%  provided by objects derived from the Function class.
%    Note that if the state signal has the reference filter enabled, the terminal
%  penalty becomes || Q (x_N - x_ref) ||_p, where x_ref  is a user-specified
%  reference signal, taken from model.x.reference.
%    To add this filter, call model.x.with('terminalPenalty'). Then you can specify
%  parameters of the penalty function by setting the model.x.terminalPenalty
%  property to an instance of the Function class (see " help Function").
%    To remove this filter, call model.x.without('terminalPenalty').
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
 
 
if ~obj.isKind('x')
	error('Filter "terminalPenalty" can only be added to state variables.');
end

% set up the filter
filter = FilterSetup;
filter.addField('value', []);

% the filter impacts the following calls:
filter.callback('objective') = @on_objective;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when constructing the objective function

out = 0;
if isempty(obj.terminalPenalty)
	return
end

% no reference by default
reference = zeros(obj.n, 1);

if obj.hasFilter('reference')
	% reference can either be free (sdpvar) or fixed (last column)
	if ismember(obj.Internal.reference.type, {'free', 'symbolic'})
		% symbolic reference, we implicitly assume it's a vector
		reference = obj.Internal.reference.var;
	elseif ~isempty(obj.reference)
		% fixed reference
		reference = obj.reference(:, end);
	end
end

out = obj.terminalPenalty.feval(obj.var(:, end) - reference);

end

%------------------------------------------------
function obj = on_set(obj, P)
% called prior to property being set

% validate the penalty (empty penalty means no penalization)
if ~isempty(P)
	error(obj.validatePenalty(P));
end
obj.terminalPenalty = P;

end
