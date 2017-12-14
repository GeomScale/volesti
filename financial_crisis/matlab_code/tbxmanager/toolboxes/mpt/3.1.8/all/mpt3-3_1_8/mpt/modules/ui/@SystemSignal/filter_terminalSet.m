function filter = filter_terminalSet(obj)
%
%  TERMINALSET: Adds a polyhedral constraint on the final predicted state 
%  =======================================================================
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
%     This filter adds a polyhedral constraint on the final predicted state, i.e.,
%  x_N in P, where x_Nis the final predicted state and P  is a polyhedron of
%  suitable dimension.
%    To enable the filter, first call model.x.with('terminalSet'). Then you can
%  provide the polyhedron in the model.x.terminalSet property.
%    To remove this filter, call model.x.without('terminalSet').
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
	error('Filter "terminalSet" can only be added to state variables.');
end

% set up the filter
filter = FilterSetup;
filter.addField('value', [], @validate_polyhedron);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

% exit immediately if no terminal set is provided
if ~isa(obj.terminalSet, 'Polyhedron')
	out = [];
	return
end

H = obj.terminalSet.A;
K = obj.terminalSet.b;
Heq = obj.terminalSet.Ae; 
Keq = obj.terminalSet.be;

% TODO: support multiple terminal sets
out = [ H*obj.var(:, end) <= K ];
if ~isempty(Heq)
    out = out + [ Heq*obj.var(:, end) == Keq ];
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called before the terminal set is changed

% empty penalty means no penalization
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
Q = P.copy(); % make a copy
for i = 1:numel(Q)
	Q(i).minHRep().normalize();
end

obj.terminalSet = Q;

end
