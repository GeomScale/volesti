function filter = filter_min(obj)
%
%  MIN: Lower bound on a signal 
%  =============================
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
%     This filter introduces the lower bound lb  on the signal s, i.e. lb <= s_k .
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
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

out = [];
for i = 1:obj.n
	% Do not include +/-Inf bounds
	if any(~isinf(obj.min(i, :)))
		out = out + [ obj.min(i, :) <= obj.var(i, :) ];
	end
end

end

%------------------------------------------------
function obj = on_set(obj, value)
% called when the filter's values are changed

value = value(:);
if numel(value)~=obj.n
	error('Value must be a %dx1 vector.', obj.n);
end
if obj.hasFilter('max') && any(value > obj.max)
	error('Lower bound cannot exceed upper bound.');
end
obj.min = value;

end
