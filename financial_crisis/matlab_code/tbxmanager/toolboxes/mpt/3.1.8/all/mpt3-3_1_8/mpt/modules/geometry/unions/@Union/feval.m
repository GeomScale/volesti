function [fval, feasible, idx, tie_value] = feval(obj, x, varargin)
%
%  FEVAL: Evaluates a given function defined over a union of convex sets. 
%  =======================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      fval = U.feval(x)
%      fval = U.feval(x, function_name)
%      [fval, feasible, idx, tb_value] = U.feval(x, function_name)
%      [fval, feasible, idx, tb_value] = feval(U, x, function_name)
%    
%  
%  DESCRIPTION
%  -----------
%     Evaluates function for a given value of the point x over the union of convex
%  sets U characterized by the name function_name. If the string function_name is
%  omitted, it is assumed that only one function is attached to the union. The
%  dimension of the vector x must be the same as the dimension of all sets in the
%  union. The evaluation is based on the following approach: 
%    
%     1. Find the index of the set where the point x lies. 
%     2. Evaluate the function over the set determined by given index. 
%    If the point lies outside of the union, the result is NaN. Notes:
%     
%    
%     - U must be a single union. Arrays of unions are not accepted. Use
%     array.forEach(@(e) e.feval(...)) to evaluate arrays of unions. 
%     - 
%   function_name must refer to a single function. If omitted, U.feval(x) only
%     works if the union has a single function. 
%     Outputs
%     
%    
%     1. x is not contained in any set of the union:
%   fval = (m X 1) vector of NaNs, where m is the range of the function, feasible =
%     false, idx = [], tb_value = []. 
%     2. 
%   x is in a single set:
%   fval = (m X 1) vector of function values, feasible = true, idx = index of the
%     set which contains x, tb_value = []. 
%     3. 
%   x is contained in multiple sets (either at the boundary or in strict interior
%     if there are overlaps), no tie-breaking (default):
%   fval = (m X j) matrix of function values (j denotes the number of sets which
%     contain x), each column contains the value of function_name in the
%     corresponding set, feasible = true, idx = (1 X j) vector of indices of sets
%     which contain x, tb_value = []. 
%     4. x is contained in multiple sets (either at the boundary or in strict
%     interior if there are overlaps), tie-breaking enabled (see below):
%   fval = (m X 1) vector containing the function value in the set in which value
%     of the tie-breaking function is smallest (if there are multiple sets with the
%     same tie-breaking value, the first such set is considered), feasible = true,
%     idx = index of the set which contains x and, simultaneously, has the smallest
%     value of the tie-breaking function, tb_value = scalar value of the
%     tie-breaking function in set indexed by idx. 
%     Tie-breaking
%    The purpose of tie-breaking is to automatically resolve situations where the
%  evaluation point x is contained in multiple sets. With tie-breaking enabled
%  Union/feval() evaluates the tie-breaking function to decide which set containing
%  x should be used for evaluation of the original function. The tie-breaking
%  function can be specified by U.feval(x, 'tiebreak', tb_fun), where tb_fun can be
%  either a string or a function handle. A string value must refer to an another
%  function which exists in the union U. A typical case where tie-breaking is
%  useful is evaluation of discontinuous MPC feedback laws: uopt = U.feval(x,
%  'primal', 'tiebreak', 'obj')Here, if x is contained in multiple sets, then the
%  function primalis only evaluated in the set which contain x and simultaneously
%  has the smallest value of the tie-breaking function obj. A special case of
%  tie-breaking is the "first-set" rule where we are only interested in evaluating
%  a function in the first set which contains x (despite the fact there can be
%  multiple such sets). This is achieved by fval = U.feval(x, 'function_name',
%  'tiebreak', @(x) 0)Notes: 
%    
%     - Tie-breaking functions must be scalar-valued. 
%     - No tie-breaking is done by default. 
%     Evaluation in particular sets
%     fval = U.feval(x, 'myfun', 'regions', indices) evaluates function myfun over
%  all sets indexed by indices. The output fval is always an (m X j)  matrix, where
%  j  is the cardinality of indices. Note that once the regions option is enabled,
%  Union/feval() will not perform point location. Instead, it will evaluate the
%  function in all sets indexed by indices, regardless of whether they contain x or
%  not. The regions option allows to quickly evaluate multiple functions as
%  follows: [first_value, idx] = U.feval(x, 'first_function') second_value =
%  U.feval(x, 'second_function', 'regions', idx)In the second call, Union/feval
%  will only evaluate second_functionin sets specified by the indices option, hence
%  skipping expensive point location.
%  
%  INPUT
%  -----
%     
%        
%          U             Union of convex sets derived from the    
%                        ConvexSet class, e.g. Polyhedron, YSet,  
%                        ...                                      
%                        Class: Union                             
%          x             A point at which the function should be  
%                        evaluated. The point must be given as a  
%                        column real vector with the same         
%                        dimension as the convex set.             
%                        Class: double                            
%          function_name Name of the function to evaluate. The    
%                        string must match one of the stored      
%                        function names. If there is only one     
%                        function attached, this argument can be  
%                        omitted.                                 
%                        Class: char                              
%                          
%  
%  
%  OUTPUT
%  ------
%     
%        
%          fval     Function value at the point x over the   
%                   union of convex sets U.                  
%                   Class: double                            
%          feasible Logical value indicating if the point x  
%                   is contained in the union or not.        
%                   Class: logical                           
%                   Allowed values:                          
%                                                            
%                     1                                      
%                     0                                      
%                                                            
%          idx      Vector of indices of sets that contain   
%                   the point x.                             
%                   Class: double                            
%          tb_value Value of the tie-breaking function if    
%                   the point belongs to multiple sets.      
%                   Class: double                            
%                     
%  
%  
%  SEE ALSO
%  --------
%     fplot
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
narginchk(2, Inf);

fval = [];
feasible = false;
idx = [];
tie_value = [];
if numel(obj)==0
	return
end

%% parse inputs
function_name = '';
options.tiebreak = []; % no tie-breaking by default
options.regions = [];
if nargin>2
	if mod(nargin, 2)==0
		% U.feval(x, ['option', value, ...])
		start_idx = 1;
	else
		% U.feval(x, 'function', ['option', value, ...])
		function_name = varargin{1};
		start_idx = 2;
	end
	for i = start_idx:2:length(varargin)
		options.(varargin{i}) = varargin{i+1};
	end
end

%% deal with arrays
if numel(obj)>1
	% evaluate the function in each partition and optionally apply
	% tiebreaking
	fval = {};
	feasible = false;
	regs = {};
	tie_value = Inf(1, numel(obj));
	for i = 1:numel(obj)
		% evaluate the function in each partition
		[fval{i}, fs, regs{i}] = obj(i).feval(x, varargin{:});
		feasible = feasible | fs;
		if fs && ~isempty(options.tiebreak)
			% evaluate the tiebreaking function in the "active" region
			if ischar(options.tiebreak)
				tie_value(i) = obj(i).Set(regs{i}(1)).Functions(options.tiebreak).feval(x);
			else
				tie_value(i) = feval(options.tiebreak, x);
			end
		end
	end
	if ~feasible
		fval = NaN(size(fval{1}, 1), 1);
		idx = [];
		tie_value = [];
	elseif ~isempty(options.tiebreak)
		% pick the partition with smallest value of the tiebreak
		[~, ipart] = min(tie_value);
		idx = [ipart*ones(1, length(regs{ipart})); regs{ipart}];
		fval = fval{ipart};
		tie_value = tie_value(ipart);
	else
		% no tiebreaking
		fv = fval;
		fval = [];
		idx = [];
		for i = 1:numel(obj)
			if any(~isnan(fv{i}))
				fval = [fval, fv{i}];
			end
			idx = [idx, [i*ones(1, length(regs{i})); regs{i}]];
		end
		tie_value = [];
	end
	return
end

if numel(obj.Set)==0
	% no regions, nothing to do
	return
end

%% validate arguments
validate_realvector(x);
if isempty(function_name)
	fnames = obj.listFunctions();
	if isempty(fnames)
		error('The object has no functions.');
	elseif length(fnames)>1
		error('The object has multiple functions, specify the one to evaluate.');
	else
		function_name = fnames{1};
	end
elseif ~ischar(function_name)
	error('The function name must be a string.');
elseif ~obj.hasFunction(function_name)
	error('No such function "%s" in the object.', function_name);
end
if ~(isempty(options.tiebreak) || isa(options.tiebreak, 'char') || ...
		isa(options.tiebreak, 'function_handle'))
	error('The tiebreak option must be a string or a function handle.');
end
if ~isnumeric(options.regions)
	error('The regions option must be a vector of integers.');
end
if (ischar(options.tiebreak) && ~obj.hasFunction(options.tiebreak))
	error('No such function "%s" in the object.', options.tiebreak);
end


%% evaluate
iscell_set = iscell(obj.Set);

	function out = eval_region(ridx, fun)
		% applies either {i} or (i) indexing of the set, evaluates function
		% 'fun' at 'x'
		%
		% since this is an inline function, "obj" and "x" are taken from
		% the main function's workspace
		if iscell_set
			out = obj.Set{ridx}.Functions(fun).feval(x);
		else
			out = obj.Set(ridx).Functions(fun).feval(x);
		end
	end

% size of the output
m = size(eval_region(1, function_name), 1);

% which regions contain "x"?
if isempty(options.regions)
	[feasible, idx] = obj.contains(x);
else
	feasible = true;
	idx = options.regions;
end
if feasible
	% x in regions indexed by idx
	n = numel(idx);
	fval = zeros(m, n);
	if n>1 && ~isempty(options.tiebreak) && isempty(options.regions)
		% tie-breaking (not applied if regions are provided)
		tval = zeros(1, n);
		for i = 1:n
			% evaluate the tie-break function in each region
			if ischar(options.tiebreak)
				% assumes "tiebreak" is a function of the union
				tb_val = eval_region(idx(i), options.tiebreak);
			else
				% tiebreak is a function handle
				tb_val = options.tiebreak(x);
			end
			if ~isscalar(tb_val)
				error('The tie breaker must be a scalar-valued function.');
			end
			tval(i) = tb_val;
		end
		% pick the region in which the tiebreak has minimal value
		[tie_value, tie_region] = min(tval);
		idx = idx(tie_region);
		% return single function value if we have tiebreak rules
		fval = eval_region(idx, function_name);
	else
		% no tie-breaking, just evaluate the function in all regions which
		% contain "x"
		for i = 1:n
			fval(:, i) = eval_region(idx(i), function_name);
		end
	end
else
	% not in the domain
	fval = NaN(m, 1);
end

end
