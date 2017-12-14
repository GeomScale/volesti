function filter = filter_penalty(varargin)
%
%  PENALTY: Penalizes the signal in the MPC cost function 
%  =======================================================
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
%     This filter, which is enabled by default, allows to penalize a particular
%  signal in the MPC cost function. The penalty function can be specified by
%  setting the model.signal.penaltyproperty to an instance of the Function class
%  (see " help Function" for more information). Penalization of state variables is
%  achieved by setting model.x.penalty, while model.u.penalty and
%  model.y.penaltyspecify penalization of input and output variables, respectively.
%  
%    Note that this filter only adds penalties to signals predicted at steps k = 0,
%  ..., N-1  of the prediction horizon. Therefore setting model.x.penalty will NOT
%  penalize the final predicted state (i.e., x_N). To add a terminal state penalty,
%  use the terminalPenalty filter (see " help
%  SystemSignal/filter_terminalPenalty").
%    To disable penalization of a signal, set its penaltyproperty to an empty
%  matrix.
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

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when constructing the cost function

out = 0;
if isempty(obj.penalty)
	return
end

% remember that state variables have length N+1, but we only
% penalize the first N components. the terminal penalty can be
% easily added by obj.with('terminalPenalty')
if obj.isKind('x')
	M = obj.N-1;
else
	M = obj.N;
end

reference = zeros(obj.n, M);
if obj.hasFilter('reference')
	if ismember(obj.Internal.reference.type, {'free', 'preview', 'symbolic'})
		% symbolic reference (vector or matrix)
		reference = obj.Internal.reference.var;
	elseif ~isempty(obj.reference)
		% numerical reference (vector or matrix)
		reference = obj.reference;
	end
end

for k = 1:M
	% if "k" exceeds number of references, repeat with the last provided
	% reference. If size(reference,2)>1, then we get trajectory preview.
	k_ref = min(k, size(reference, 2));
	value = obj.var(:, k) - reference(:, k_ref);
	out = out + obj.penalty.feval(value);
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called prior to property being set

% validate the penalty (empty penalty means no penalization)
if ~isempty(P)
	error(obj.validatePenalty(P));
end
obj.penalty = P;

end
