function filter = filter_binary(obj)
%
%  BINARY: Constraints variable to be binary (0/1) 
%  ================================================
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
%     Adding this filter will constraint some (or all) elements of a given variable
%  to take only binary values. To enable the filter, use
%  model.signal.with('binary').
%    To impose binary on all elements of a given variable (say, model.u), use
%  model.u.binary = true. To add binary only to elements indexed by idx, call
%  model.u.binary = idx. To mark all elements of model.u as real variables, use
%  model.u.binary = [].
%    To remove this filter, call model.signal.without('binary'), in which case all
%  elements of signal will be considered as real-valued variables.
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
filter.addField('index', 1:obj.n, @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('set') = @on_set;

end


%------------------------------------------------
function obj = on_set(obj, new)
% called when the filter's parameters are changes

if islogical(new) && new
	% true = all elements are binary
	obj.binary = 1:obj.n;
elseif islogical(new) && ~new
	% false = no binary elements
	obj.binary = [];
else
	% just set indices
	obj.binary = new;
end

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

if islogical(obj.binary)
	idx = 1:obj.n;
else
	idx = obj.binary;
end

out = [ binary(obj.var(obj.binary, :)) ];
out = out + [ 0 <= obj.var(obj.binary, :) <= 1 ];

end
