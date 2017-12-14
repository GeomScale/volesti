classdef ConvexSet < ConvexSetInterface & IterableBehavior & matlab.mixin.Copyable
%
%  CONVEXSET: Represets a convex set in MPT 
%  =========================================
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
%     Represents a general convex set in R^n. The dimension of the set, i.e. nis
%  stored under Dim property. The ConvexSet cannot be instantiated, you can only
%  create objects derived from this class, see related functions. You can associate
%  functions to ConvexSet class via addFunction methods. Function handles are
%  stored inside Func property. The ConvexSet class handles various methods that
%  operate over convex sets and functions over convex sets. For a list of available
%  methods type "methods('ConvexSet')".
%  
%  SEE ALSO
%  --------
%     YSet,  Polyhedron,  Union,  PolyUnion
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
  properties (SetAccess=protected)
      Dim % Dimension of the set
  end
  properties (SetAccess=protected, Hidden)
      Internal % Internal data
	  Functions % Storage of attached functions
  end

  properties
      Data % any user changeable data
  end
  
  properties (Dependent=true, Transient=true, Hidden)
	  Func
  end
  
  methods
	  
	  function F = fliplr(P)
		  % Flips an array of ConvexSet objects

		  % since we store arrays as column vectors by default (see
		  % ConvexSet/horzcat), it suffices to invoke flipud()
		  F = flipud(P);
	  end

	  function F = get.Func(obj)
		  F = obj.Functions.values();
	  end
	  
	  function obj = addFunction(obj, F, FuncName)
		  % adds a function to the object
		  
		  narginchk(3, 3);
		  if numel(obj)==1
			  % check domain
			  if (isa(F, 'AffFunction') || isa(F, 'QuadFunction')) && ...
					  obj.Dim ~= F.D
				  error('The function must have the same domain as the set.');
			  elseif isa(F, 'function_handle')
				  F = Function(F);
			  elseif isa(F, 'char')
				  error('First input must be a function object.');
			  end
			  % add the function to the map
			  obj.Functions(FuncName) = F;
			  
		  elseif numel(obj)>1
			  % deal with arrays
			  for j=1:numel(obj)
				  obj(j) = obj(j).addFunction(F,FuncName);
			  end
		  end
	  end
	  
	  function F = getFunction(obj, FuncName)
		  % returns function indexed by the string FuncName
          if ~isa(FuncName,'char')
              error('The function name must be given as a string.');
          end
          for i = 1:numel(obj)
              F(i) = obj(i).Functions(FuncName);
          end
	  end
	  
	  function obj = removeFunction(obj, FuncNames)
		  % removes function indexed by the string FuncName
          if iscell(FuncNames)
               if any(~cellfun('isclass',FuncNames,'char'))
                   error('The function name must be given as a string.');
               end
          elseif ~isa(FuncNames,'char')
              error('The function name must be given as a string.');
          end
		  for i = 1:numel(obj)
			  obj(i).Functions.remove(FuncNames);
		  end
	  end
	  
	  function obj = removeAllFunctions(obj)
		  % removes all attached functions
		  for i = 1:numel(obj)
			  obj(i).Functions.remove(obj(i).Functions.keys);
		  end
	  end
	  
	  function out = listFunctions(obj)
		  % lists attached functions
		  %
		  % outputs:
		  % * empty cell array if "obj" is an empty object
		  % * cell array of function names if "obj" is a single set
		  % * cell array of cell arrays of function names if "obj" is an
		  %   array

		  if numel(obj)==0
			  out = {};
		  elseif numel(obj)==1
			  out = obj.Functions.keys();
		  else
			  out = cell(1, numel(obj));
			  for i = 1:numel(obj)
				  out{i} = obj(i).Functions.keys();
			  end
		  end
	  end
	  
	  function out = hasFunction(obj, FuncName)
		  % returns true if the object contains function(s) indexed by
		  % FuncName
		  %
		  % inputs:
		  %   FuncName: either a string or a cell array of strings
		  %
		  % outputs:
		  % * empty double if "obj" is an empty object
		  % * column logical vector if "obj" is a single set (each row
		  %   corresponds to presence of FuncName{i})
		  % * logical matrix if "obj" is an array, with "n" rows and "m"
		  %   columns, where "n" is the number of functions which are
		  %   queried, and "m" is the number of sets. then out(i, j)=true
		  %   if the j-th set contains FuncName{i}

		  if numel(obj)==0
			  out = [];
		  elseif numel(obj)==1
			  x = obj.Functions.isKey(FuncName);
			  % make sure we always return column vector if we have
			  % multiple functions
			  out = x(:); 
		  else
			  if iscell(FuncName)
				  n = length(FuncName);
			  else
				  n = 1;
			  end
			  out = false(n, numel(obj));
			  for i = 1:numel(obj)
				  x = obj(i).Functions.isKey(FuncName);
				  out(:, i) = x(:);
			  end
		  end
	  end

	  function [F, map] = uniqueFunctions(obj, FuncName)
		  % Returns unique occurences of function "FuncName" and their map
		  % to elements of the array.
		  %
		  % Given an array of sets "P" whose all elements have the function
		  % "FuncName" attached, this method returns list of unique
		  % functions "F" attached to each element of "P". The method also
		  % returns indices "map", such that P(i) uses function F(map(i)).

		  F = []; % list of unique functions
		  map = zeros(1, numel(obj));
		  for i = 1:numel(obj)
			  if ~obj(i).hasFunction(FuncName)
				  error('No such function "%s".', FuncName);
			  end
			  f = obj(i).getFunction(FuncName);
			  % is the function unique?
			  u = find(arrayfun(@(x) x==f, F));
			  if isempty(u)
				  % unique function found
				  F = [F f];
				  map(i) = length(F);
			  else
				  % another identical function was already seen
				  map(i) = u(1);
			  end
		  end
	  end

	  function out = trimFunction(obj, FuncName, n)
		  % Extracts the first "n" rows of a given affine function
		  %
		  % This method creates a copy of the input array where the
		  % specified affine function is replaced by a new affine function
		  % which only contains the first "n" rows of the original
		  % function.
		  
		  narginchk(3, 3);
		  [~, errmsg] = obj.validateFunctionName(FuncName);
		  error(errmsg);
		  
		  if nargout==0
			  % in-place trimming
			  out = obj;
		  else
			  % create a copy if explicitly requested
			  out = obj.copy();
		  end

		  for i = 1:numel(out)
			  f = out(i).getFunction(FuncName);
			  if ~isa(f, 'AffFunction')
				  error('Only affine functions can be trimmed.');
			  end
			  ft = AffFunction(f.F(1:n, :), f.g(1:n), f.Data);
			  out(i).addFunction(ft, FuncName);
		  end
	  end
	  
  end
  
  methods (Hidden)
	  % private APIs
  
	  function target = copyFunctionsFrom(target, source)
		  % copies functions attached to "source" to "target"
		  
		  if length(source.Functions)>0
			  funs = source.listFunctions;
			  for i=1:length(funs)
				  target.addFunction(source.getFunction(funs{i}), funs{i});
			  end
		  end

	  end
	  
      function obj = setInternal(obj,name,value)
          %
          % overwrites internal property for the ConvexSet object
          % (internal function)
          %
          % If we want to add the internal property from outside of
          % this class (e.g. inside the PLCP solver) e.g.
          %
          % obj.Internal.name = value,
          %
          % use the syntax:
          %         obj.setInternal('name',value)
          %
          % DO NOT USE THIS METHOD UNLESS YOU PERFECTLY KNOW WHAT YOU ARE DOING          
          %
          
          narginchk(3, 3);
          
          if ~ischar(name)
              error('Name must be a string.');
          end
          
          obj.Internal.(name) = value;
          
      end
  end
  
  methods (Sealed)
	  
	  % plotting dispatchers (actual plotting is done by plot_internal and
	  % fplot_internal)
	  h = fplot(obj, varargin)
	  h = plot(varargin)
	  
  end
  
  methods (Access=protected)

	  % function prototypes (plotting methods must remain protected)
	  h = fplot_internal(obj, function_name, options)
	  h = plot_internal(obj, options)

	  function new = copyElement(obj)
		  % Create a copy of an object

		  % shallow copy of obj.Internal and obj.Data
		  new = copyElement@matlab.mixin.Copyable(obj);
		  
		  % deep copy of obj.Functions
		  keys = obj.Functions.keys;
		  if isempty(keys)
			  new.Functions = containers.Map;
		  else
			  new.Functions = containers.Map(keys, obj.Functions.values);
		  end
		  		  
	  end
	  
	  function displayFunctions(obj)
		  % displays attached functions (to be used from a disp() method)
		  
		  fprintf('Functions');
		  funs = obj.listFunctions;
		  nf = length(funs);
		  if nf==0
			  fprintf(' : none\n');
		  else
			  fprintf(' : %d attached ',nf);
			  for i = 1:nf
				  fprintf('"%s"', funs{i});
				  if i<nf
					  fprintf(', ');
				  end
			  end
			  fprintf('\n');
		  end
	  end

	  function [function_name, msg] = validateFunctionName(obj, function_name)
		  % returns a valid function_name if none was provided
		  % returns an error message if validation fails

		  msg = '';
		  if isempty(function_name)
			  fnames = obj(1).listFunctions();
			  if isempty(fnames)
				  msg = 'The object has no functions.';
				  return
			  elseif numel(fnames)>1
				  msg = 'The object has multiple functions, specify the one to evaluate.';
				  return
			  else
				  function_name = fnames{1};
			  end
			  % check that all remaining sets have this function defined
			  for i = 2:numel(obj)
				  if ~obj(i).hasFunction(function_name)
					  msg = sprintf('No such function "%s" in set %d.', function_name, i);
					  return
				  end
			  end
		  elseif ~ischar(function_name)
			  msg = 'The function name must be a string.';
			  return
		  elseif any(~obj.hasFunction(function_name))
			  msg = sprintf('No such function "%s" in the object.', function_name);
			  return
		  end

	  end

  end
  
end

