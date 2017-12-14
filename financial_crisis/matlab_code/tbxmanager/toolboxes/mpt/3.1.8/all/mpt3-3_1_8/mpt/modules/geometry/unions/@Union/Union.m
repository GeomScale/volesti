classdef Union < handle & IterableBehavior & matlab.mixin.Copyable
%
%  UNION: Represents a general union of convex sets in MPT 
%  ========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U = Union(Set)
%    
%  
%  DESCRIPTION
%  -----------
%     The Union object represent collection of various convex sets. The only
%  restriction for the sets is to be convex, i.e. they have to be derived from the
%  ConvexSet class. You can associate functions to any of the set via addFunction
%  method of the ConvexSet class. Function handles and all properties of the sets
%  can be accessed via Union.Set property based on the index. For a list of
%  available methods type "methods('Union')".
%  
%  INPUT
%  -----
%     
%        
%          Set Any object derived from the ConvexSet    
%              class.                                   
%              Class: ConvexSet                         
%                
%  
%  
%  OUTPUT
%  ------
%     
%        
%          U Object of the Union class.               
%            Class: Union                             
%              
%  
%  
%  SEE ALSO
%  --------
%     YSet,  Polyhedron,  PolyUnion
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
 
 
  properties (SetAccess=private, Dependent=true)
      Num % number of sets
  end
  properties (SetAccess=protected)
      Set % container for sets
	  Domain % domain of the union
  end
  properties (SetAccess=protected, Hidden)
      Internal % internal data
  end
  properties
      Data % any user changeable data
  end
  
  methods
      function obj = Union(varargin)
          % 
          % constructor of Union object
          %
          
          % short syntax
          if nargin==1
              arg{1} = 'Set';
              arg{2} = varargin{1};
          else
              arg = varargin;
          end
          
          % full syntax
          ip = inputParser;
          ip.KeepUnmatched = false;
          ip.addParamValue('Set', [], @(x) isa(x, 'ConvexSet'));
		  ip.addParamValue('Domain', [], @(x) isa(x, 'ConvexSet'));
          ip.addParamValue('Data', [], @(x) true);
          ip.parse(arg{:});
          p = ip.Results;

          % remove empty sets
          if ~builtin('isempty',p.Set)
              C = p.Set(:);
              c = isEmptySet(C);
              C(c) = [];
          else
              C = p.Set;
          end
          
          % check attached functions, if they are the same in all sets
		  if numel(C)>0
			  funs = C(1).listFunctions();
			  for i = 2:numel(C)
				  if any(~C(i).hasFunction(funs))
					  error('All sets must have associated the same number of functions.');
				  end
			  end
		  end

          % assign data
          obj.Set = num2cell(C(:));
          obj.Data = p.Data;
		  if ~isempty(p.Domain) && ~isempty(C)
			  % set the domain if the set is not empty
			  obj.Domain = num2cell(p.Domain(:));
		  end
	  end
	  
	  function D = get.Domain(obj)
		  % getter for Union.Domain
		  
		  if isempty(obj.Domain)
			  D = obj.Set;
		  else
			  D = obj.Domain;
		  end
	  end
	 
	  function obj = addFunction(obj, fun, FuncName)
		  % adds a function to each member of the union
		  
		  narginchk(3, 3);
		  
		  % attach the function to each element of the set
		  for i = 1:numel(obj)
			  cell_set = iscell(obj(i).Set);
			  for j = 1:length(obj(i).Set)
				  if cell_set
					  obj(i).Set{j}.addFunction(fun, FuncName);
				  else
					  obj(i).Set(j).addFunction(fun, FuncName);
				  end
			  end
		  end
	  end
		  
		  
	  function U = getFunction(obj, FuncName)
		  % returns function indexed by the string FuncName

		  narginchk(2, 2);
		  for i = 1:numel(obj)
			  % make sure the function exists
			  if any(~obj(i).hasFunction(FuncName))
				  if iscell(FuncName)
					  idx = find(~obj(i).hasFunction(FuncName));
					  missing = FuncName{idx(1)};
				  else
					  missing = FuncName;
				  end
				  error('No such function "%s" in the object.', missing);
			  end
			  
			  U(i) = obj(i).copy();
			  toremove = setdiff(obj(i).listFunctions, FuncName);
			  if ~isempty(toremove)
				  U(i).removeFunction(toremove);
			  end
		  end
	  end
	  
	  function obj = removeFunction(obj, FuncNames)
		  % removes function indexed by the string FuncName

		  narginchk(2, 2);
		  % make a copy before removing function(s)
		  for i = 1:numel(obj)
			  % make sure the function exists
			  if any(~obj(i).hasFunction(FuncNames))
				  if iscell(FuncNames)
					  idx = find(~obj(i).hasFunction(FuncNames));
					  missing = FuncNames{idx(1)};
				  else
					  missing = FuncNames;
				  end
				  error('No such function "%s" in the object.', missing);
			  end
			  if iscell(obj(i).Set)
				  for j = 1:length(obj(i).Set)
					  obj(i).Set{j}.removeFunction(FuncNames);
				  end
			  else
				  obj(i).Set.removeFunction(FuncNames);
			  end
		  end
	  end
	  
	  function obj = removeAllFunctions(obj)
		  % removes all attached functions

		  for i = 1:numel(obj)
			  if iscell(obj(i).Set)
				  for j = 1:length(obj(i).Set)
					  obj(i).Set{j}.removeAllFunctions;
				  end
			  else
				  obj(i).Set.removeAllFunctions;
			  end
		  end
	  end
	  
	  function out = listFunctions(obj)
		  % lists attached functions
		  %
		  % outputs:
		  % * empty cell array if "obj" is an empty object
		  % * cell array of function names if "obj" is a single union
		  % * cell array of cell arrays of function names if "obj" is an
		  %   array

		  if numel(obj)==0 || numel(obj(1).Set)==0
			  out = {};
		  elseif numel(obj)==1
			  if iscell(obj.Set)
				  out = obj.Set{1}.Functions.keys();
			  else
				  out = obj.Set(1).Functions.keys();
			  end
		  else
			  out = cell(1, numel(obj));
			  for i = 1:numel(obj)
				  if iscell(obj(i).Set)
					  out{i} = obj(i).Set{1}.Functions.keys();
				  else
					  out{i} = obj(i).Set(1).Functions.keys();
				  end
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
		  % * column logical vector if "obj" is a single union (each row
		  %   corresponds to presence of FuncName{i})
		  % * logical matrix if "obj" is an array, with "n" rows and "m"
		  %   columns, where "n" is the number of functions which are
		  %   queried, and "m" is the number of unions. then out(i, j)=true
		  %   if the j-th union contains FuncName{i}

		  narginchk(2, 2);
		  out = []; % default output for empty arrays
		  if numel(obj)==1 && numel(obj.Set)>0
			  if iscell(obj.Set)
				  x = obj.Set{1}.Functions.isKey(FuncName);
			  else
				  x = obj.Set(1).Functions.isKey(FuncName);
			  end
			  % make sure we always return column vector if we have
			  % multiple functions
			  out = x(:);
		  elseif numel(obj)>1
			  if iscell(FuncName)
				  n = length(FuncName);
			  else
				  n = 1;
			  end
			  out = false(n, numel(obj));
			  for i = 1:numel(obj)
				  if numel(obj(i).Set)>0
					  if iscell(obj(i).Set)
						  x = obj(i).Set{1}.Functions.isKey(FuncName);
					  else
						  x = obj(i).Set(1).Functions.isKey(FuncName);
					  end
					  out(:, i) = x(:);
				  end
			  end
		  end
	  end

	  function U = trimFunction(obj, FuncName, n)
		  % Extracts the first "n" rows of a given affine function
		  %
		  % This method creates a copy of the input union where the
		  % specified affine function is replaced by a new affine function
		  % which only contains the first "n" rows of the original
		  % function.
		  
		  if length(obj)>1
			  error('Single union please.');
		  elseif ~obj.hasFunction(FuncName)
			  error('No such function "%s".', FuncName);
		  end
		  
		  if nargout==0
			  % in-place trimming
			  U = obj;
		  else
			  % create a copy if explicitly requested
			  U = obj.copy;
		  end

		  for i = 1:U.Num
			  f = U.index_Set(i).getFunction(FuncName);
			  if ~isa(f, 'AffFunction')
				  error('Only affine functions can be trimmed.');
			  end
			  ft = AffFunction(f.F(1:n, :), f.g(1:n), f.Data);
			  U.Set(i).addFunction(ft, FuncName);
		  end
		  
	  end
	  
	  function map = findUnique(obj, function_name, varargin)
		  % Finds regions of the underlying set which share the same
		  % expression of a given function
		  %
		  % Syntax:
		  %   map = obj.findUnique(fun)
		  %   map = obj.findUnique(fun, 'range', R)
		  %
		  % Inputs:
		  %   obj: Union object
		  %   fun: String representation of the function to use
		  %     R: When specified, only the R components of a given
		  %        function are employed to recognize unique regions
		  %
		  % Output:
		  %   map: 1xN array (N=number of regions in the union's underlying
		  %        set) with map(i)=j if the i-th region contains the j-th
		  %        unique expression of function FUN
		  
		  %% validation
		  error(obj.rejectArray());
		  narginchk(2, Inf);
		  if ~ischar(function_name)
			  error('The function name must be a string.');
		  elseif ~obj.hasFunction(function_name)
			  error('No such function "%s" in the object.', function_name);
		  elseif ~isa(obj.index_Set(1).Functions(function_name), 'AffFunction')
			  error('Function "%s" must be affine.', function_name);
		  end
		  
		  %% parsing
		  ip = inputParser;
		  ip.addParamValue('range', 1:obj.index_Set(1).Functions(function_name).R, ...
			  @validate_indexset);
		  ip.parse(varargin{:});
		  options = ip.Results;

		  % clone the union since trimming changes the original
		  new = obj.copy();
		  new.trimFunction(function_name, options.range);
		  [~, map] = new.Set.uniqueFunctions(function_name);
	  end
	  
  end
  methods (Hidden)

      function obj = setInternal(obj,name,value)
          %
          % overwrites internal property for the Union object
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
  
  methods (Access=protected)
	  % internal APIs

	  % function prototypes (plotting methods must remain protected)
	  h = plot_internal(obj, idx, varargin)

	  function U = copyElement(obj)
		  % Creates a copy of the union
		  %
		  %   copy = U.copy()

		  % Note: matlab.mixin.Copyable.copy() automatically properly
		  % copies arrays and empty arrays, no need to do it here.
		  % Moreover, it automatically creates the copy of proper class.
		  U = copyElement@matlab.mixin.Copyable(obj);

		  % deep copy of obj.Set
		  if iscell(obj.Set)
			  % per-element copying
			  U.Set = cell(size(obj.Set));
			  for i = 1:numel(obj.Set)
				  U.Set{i} = obj.Set{i}.copy();
			  end
		  else
			  % resort to ConvexSet/copy
			  U.Set = obj.Set.copy();
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
	  
	  function out = index_Set(obj, i)
		  % returns the i-th element of obj.Set, regardless of whether Set
		  % is a cell or an ordinary array

		  if iscell(obj.Set)
			  out = obj.Set{i};
		  else
			  out = obj.Set(i);
		  end
	  end

  end

  methods (Sealed)
	  
	  % plotting dispatchers (actual plotting is done by plot_internal)
	  h = plot(varargin)
	  
  end

  methods      
   function n = get.Num(obj)
       % get method for number of sets
       n = numel(obj.Set);
   end
  end
  
end

