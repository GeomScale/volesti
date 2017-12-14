function filter = filter_setConstraint(obj)
%
%  SETCONSTRAINT: Adds a polyhedral constraint 
%  ============================================
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
%     This filter adds polyhedral constraints z_k in P  for all k = 0, ..., N-1  to
%  the MPC setup. The filter can be applied to arbitrary signals, be it state,
%  input, or output variables.
%    Note, however, that adding this filter to a state variable will NOT add the
%  constraint on the final predicted state (only quantities on prediction steps k =
%  0, ..., N-1  are affected by this filter). To add a polyhedral terminal state
%  constraint, use the terminalSet filter (see " help
%  SystemSignal/filter_terminalSet").
%    To enable the filter, first use model.x.with('setConstraint') (you can replace
%  the x signal by any other signal of the prediction model).
%    The polyhedron P  of the new constraint x_k in P  can then be specified in the
%  setConstraint property of the signal:
%     model.x.setConstraint = Polyhedron(...)
%    To remove this filter, call model.x.without('setConstraint').
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
filter.addField('value', [], @validate_polyhedron);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

P = obj.setConstraint;
H = P.A;
K = P.b;
Heq = P.Ae;
Keq = P.be;

out = [];
for k = 1:obj.N
    out = out + [ H*obj.var(:, k) <= K ];
    if ~isempty(Heq)
        out = out + [ Heq*obj.var(:, k) == Keq ];
    end
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called before the set is changed

if isempty(P)
    return
elseif ~isa(P, 'Polyhedron')
	error('The input must be a polyhedron.');
elseif numel(P)~=1
	error('The input must be a single polyhedron.');
elseif P.Dim ~= obj.n
	error('The polyhedron must be in dimension %d.', obj.n);
elseif P.isEmptySet()
	error('The polyhedron must not be empty.');
end

% we require the H-representation, which we also normalize to avoid
% numerical problems
Q = Polyhedron(P); % make a coppy
for i = 1:numel(Q)
	Q(i).minHRep().normalize();
end

obj.setConstraint = Q;

end
