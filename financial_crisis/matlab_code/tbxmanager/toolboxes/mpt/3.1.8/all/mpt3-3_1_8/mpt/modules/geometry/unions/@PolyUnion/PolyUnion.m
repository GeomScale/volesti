classdef PolyUnion < Union
%
%  POLYUNION: Represents a union of polyhedra in the same dimension 
%  =================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      PU = PolyUnion(P)
%      PU = PolyUnion('Set',P,'Domain',D,'Convex',true,'Overlaps',false,...)
%    
%  
%  DESCRIPTION
%  -----------
%     The PolyUnion object represent collection of polyhedra in the same dimension.
%  The only restriction to construct the PolyUnion object is that the sets have a
%  common dimension. The union thus can contain bounded, unbounded, and
%  lower-dimensional polyhedra. Empty sets are removed from at the time of
%  construction of PolyUnion object. You can associate various properties with the
%  PolyUnion object: 
%    
%     - Domain of the union can be set via the Domain property. If not provided,
%     			the domain will be equal to the underlying set of the union. 
%     - To mark the union as convex, you can set the Convex property true or false.
%     If you specify that the union is convex, the remaining algorithms will rely
%     on this information which will can  speed up the computations significantly. 
%     - To mark that there are overlaps in the union, you can set the Overlaps
%     property as true or false. If the union contains non-overlapping polyhedra,
%     the union forms a partition. 
%     - To mark that the union is connected in the union, you can set the Connected
%     property as true or false. If the union is convex, if implies that the union
%     is also connected. 
%     - To mark that the union comprises only of bounded polyhedra, you can set the
%     Bounded property as true or false.  
%     - To mark that the union contains only full-dimensional polyhedra, you can
%     set the FullDim property as true or false.  
%    By specifying any of the above properties you promise that these properties
%  are automatically satisfied by the provided set, otherwise unexpected results
%  may occur. These properties can be accessed via PU.Internal field. You can
%  associate functions to any of the set via addFunction method of the ConvexSet
%  class. Function handles and all properties of the sets can be accessed via
%  PolyUnion.Set property based on the index. For a list of available methods type
%  "methods('PolyUnion')".
%  
%  INPUT
%  -----
%     
%        
%          P         Non-empty polyhedra in the same          
%                    dimension.                               
%                    Class: Polyhedron                        
%          Domain    Specifies domains of the union.          
%                    Class: Polyhedron                        
%          Convex    Set this property to true if you are     
%                    sure that the union of the polyhedra is  
%                    convex.                                  
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      true                                   
%                      false                                  
%                                                             
%          Overlaps  Set this property to false if you are    
%                    sure that the union of the polyhedra is  
%                    non-overlapping.                         
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      true                                   
%                      false                                  
%                                                             
%          Connected Set this property to true if you are     
%                    sure that the union of the polyhedra is  
%                    connected. If the union is convex, then  
%                    this property is automatically set as    
%                    true.                                    
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      true                                   
%                      false                                  
%                                                             
%          Bounded   Set this property to true if you are     
%                    sure that all the polyhedra in the union 
%                    are bounded.                             
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      true                                   
%                      false                                  
%                                                             
%          FullDim   Set this property to true if you are     
%                    sure that all the polyhedra in the union 
%                    are full-dimensional.                    
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      true                                   
%                      false                                  
%                                                             
%                      
%  
%  
%  OUTPUT
%  ------
%     
%        
%          PU Object of the PolyUnion class.           
%             Class: PolyUnion                         
%               
%  
%  
%  SEE ALSO
%  --------
%     Polyhedron,  Union
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
      Dim % Dimension of the union
  end
  
  methods
      function obj = PolyUnion(varargin)
          % 
          % Union(regions)
          % Union('Set',regions,'convexity',true,'overlaps',false);
          % Union('Set',regions,'convexity',true,'overlaps',false);
          
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
          ip.addParamValue('Set', [], @(x) isa(x, 'Polyhedron'));
          ip.addParamValue('Domain', [], @(x) isa(x, 'Polyhedron'));
		  ip.addParamValue('Convex',[], @(x) islogical(x) || x==1 || x==0);
          ip.addParamValue('Overlaps',[], @(x) islogical(x) || x==1 || x==0);
          ip.addParamValue('Connected',[], @(x) islogical(x) || x==1 || x==0);
          ip.addParamValue('Bounded',[], @(x) islogical(x) || x==1 || x==0);
          ip.addParamValue('FullDim',[], @(x) islogical(x) || x==1 || x==0);
          ip.addParamValue('Data', [], @(x) true);
          ip.parse(arg{:});
          p = ip.Results;
                                 
          % remove empty sets
          if ~builtin('isempty',p.Set)
              C = p.Set(:);
              c = isEmptySet(C);
              C = C(~c);
          else
              C = p.Set;
          end
                    
          nC = numel(C);
          if nC>0
              if ~isEmptySet(C)
                  % check dimension
                  D = [C.Dim];
                  if any(D(1)~=D)
                      error('All polyhedra must be in the same dimension.');
                  end
                  obj.Dim = D(1);
              end
              
              % assign properties
              obj.Set = C;
              %obj.Num = length(obj.Set);
              
              
              % check attached functions, if they are the same in all sets
			  if numel(C)>0
				  funs = C(1).listFunctions();
				  for i = 2:numel(C)
					  if any(~C(i).hasFunction(funs))
						  error('All sets must have associated the same number of functions.');
					  end
				  end
			  end
		  end

		  if ~isempty(p.Domain) && ~isempty(C)
			  % set the domain if the set is not empty
			  
			  % first kick out empty sets
			  p.Domain = p.Domain(~p.Domain.isEmptySet);
			  
			  if ~isempty(p.Domain)
				  obj.Domain = p.Domain(:);
				  if any(diff([obj.Domain.Dim]~=0))
					  error('All domains must be in the same dimension.');
				  elseif obj.Domain(1).Dim~=obj.Set(1).Dim
					  error('The domain must be a %dD polyhedron.', obj.Set(1).Dim);
				  end
			  end
		  end
		  
          obj.Internal.Convex = p.Convex;
          obj.Internal.Overlaps = p.Overlaps;
          % convex union implies connected
          if p.Convex
              obj.Internal.Connected = true;
          else
              obj.Internal.Connected = p.Connected;
          end
          obj.Internal.Bounded = p.Bounded; 
          obj.Internal.FullDim = p.FullDim; % full dimensionality
          obj.Data = p.Data;

          
	  end
	  
	  function out = findSaturated(obj, function_name, varargin)
		  % Determines saturation properties of unions
		  %
		  % Syntax:
		  %   out = obj.findSaturated(myfun, 'min', min_value, ...
		  %                          'max', max_value, 'range', range)
		  %
		  % Inputs:
		  %         obj: Union object
		  %       myfun: string name of the function to check
		  %              (the function must be affine)
		  %   min_value: lower saturation limit
		  %              (if not provided, the limit is comptued
		  %              automatically by solving an LP)
		  %   max_value: upper saturation limit
		  %              (if not provided, the limit is comptued
		  %              automatically by solving an LP)
		  %       range: range of function's outputs (optional)
		  %
		  % Output:
		  %     out.min: lower saturation limits
		  %     out.max: upper saturation limits
		  %    out.Imin: indices of regions where all RANGE elements of
		  %              MYFUN are saturated at the lower limit
		  %    out.Imax: indices of regions where all RANGE elements of
		  %              MYFUN are saturated at the upper limit
		  %  out.Iunsat: indices of unsaturated regions
		  %       out.S: saturation indicators:
		  %              out.S(i, j) = 1  if in the i-th region the j-th
		  %                     element of MYFUN attains the upper limit
		  %              out.S(i, j) = -1 if in the i-th region the j-th
		  %                     element of MYFUN attains the upper limit
		  %              out.S(i, j) = 0  if in the i-th region the j-th
		  %                     element of MYFUN is not saturated
		  
		  global MPTOPTIONS
		  
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
		  ip.addParamValue('min', [], @validate_realvector);
		  ip.addParamValue('max', [], @validate_realvector);
		  ip.parse(varargin{:});
		  options = ip.Results;
		  R = numel(options.range);
		  
		  %% find saturation limits
		  if isempty(options.min)
			  minvalue = Inf(R, 1);
			  for i = 1:obj.Num
				  fun = obj.Set(i).Functions(function_name);
				  lp = obj.Set(i).optMat;
				  lp.quicklp = true;
				  for j = 1:R
					  lp.f = fun.F(j, :);
					  sol = mpt_solve(lp);
					  if sol.exitflag ~= MPTOPTIONS.OK
						  error('Problems with solving an LP.');
					  end
					  minvalue(j) = min(minvalue(j), sol.obj+fun.g(j));
				  end
			  end
			  options.min = minvalue;
		  end
		  if isempty(options.max)
			  maxvalue = -Inf(R, 1);
			  for i = 1:obj.Num
				  fun = obj.Set(i).Functions(function_name);
				  lp = obj.Set(i).optMat;
				  lp.quicklp = true;
				  for j = 1:R
					  lp.f = -fun.F(j, :);
					  sol = mpt_solve(lp);
					  if sol.exitflag ~= MPTOPTIONS.OK
						  error('Problems with solving an LP.');
					  end
					  maxvalue(j) = max(maxvalue(j), -sol.obj+fun.g(j));
				  end
			  end
			  options.max = maxvalue;
		  end
		  
		  %% additional validation
		  error(validate_vector(options.min, R, 'min value'));
		  error(validate_vector(options.max, R, 'max value'));
		  
		  %% processing
		  S = zeros(numel(options.range), obj.Num);
		  for i = 1:obj.Num
			  fun = obj.index_Set(i).Functions(function_name);
			  for j = 1:numel(options.range)
				  jr = options.range(j);
				  % saturated function must have a zero affine term
				  if norm(fun.F(jr, :)) <= MPTOPTIONS.zero_tol
					  if abs(fun.g(jr) - options.min(j)) <= MPTOPTIONS.abs_tol
						  % region saturated at minimum
						  S(j, i) = -1;
					  elseif abs(fun.g(jr) - options.max(j)) <= MPTOPTIONS.abs_tol
						  % region saturated at maximum
						  S(j, i) = 1;
					  end
				  end
			  end
		  end
		  
		  %% create the output structure
		  out.S = S;
		  out.min = options.min;
		  out.max = options.max;
		  satmax = false(1, obj.Num);
		  satmin = false(1, obj.Num);
		  for i = 1:obj.Num
			  % which regions are saturated?
			  satmax(i) = all(S(:, i)==1);
			  satmin(i) = all(S(:, i)==-1);
		  end
		  out.Imin = find(satmin);
		  out.Imax = find(satmax);
		  out.Iunsat = setdiff(1:obj.Num, [out.Imin, out.Imax]);
	  end
	  
	  function [new, sliced_set] = slice(obj, dims, values)
		  % Slices the polyunion at given dimensions
		  %
		  % U.slice(dims, values) fixes dimensions given in "dims" to
		  % values in "values". Slices corresponding functions as well.
		  %
		  % For more information see Polyhedron/slice
		  
		  narginchk(3, 3);
          error(obj.rejectArray());
		  sliced_set = obj.Set.slice(dims, values);
		  new = PolyUnion(sliced_set);
		  % TODO: propagate properties preserved under slicing:
		  % * Convex (if original was convex, so is the slice)
		  % * Overlaps (if original didn't have them, neither will the
		  % slice)
		  % * Connected (if the original was, so will be the slice)
      end
      
      function E = envelope(obj)
          % Computes the convex envelope of a polyunion
          %
          % Given are two polyhedra P={x | A*x<=b} and Q={x | H*x<=k }. If
          % a facet of Q (and P) is not redundant for P (or Q), then it is
          % not part of the envelope.
          %
          % Technically, the envelope is the H-polyhedron E={x | At*x<=bt,
          % Ht*x<=kt} where At*x<=bt is the subsystem of A*x<=b obtained by
          % removing all the inequalities not valid for the polyhedron Q;
          % and Ht*x<=kt is defined similarly. 
          %
          % Limitations: supports only polyhedra in inequality
          % H-representation.
          %
          % Syntax:
          %   E = obj.envelope()
          %
          % Input:
          %   obj: a single PolyUnion object with an arbitrary number of
          %        polyhedra
          %
          % Output:
          %     E: convex envelope as a Polyhedron object (can contain
          %        redundant constraints. Use E.minHRep() to remove them)
          
          global MPTOPTIONS
          
          error(obj.rejectArray());
          
          % equality constraints not yet supported
          ne = obj.Set.forEach(@(x) size(x.He, 1));
          if nnz(ne)>0
              error('Equality constraints are not supported.');
          end

          % all sets must be in the minimal H-rep
          obj.Set.minHRep();

          % inequality representation of the envelope
          He = [];
          tic
          for i = 1:obj.Num
              if toc > MPTOPTIONS.report_period
                  fprintf('progress: %d/%d\n', i, obj.Num);
                  tic;
              end
              ni = length(obj.Set(i).b);
              keep = true(1, ni);
              for j = 1:obj.Num
                  if i==j, continue, end
                  for k = 1:ni
                      if ~keep(k)
                          % this constraint was already discarded
                          continue
                      end
                      % flip the k-th constraint of the i-th region and
                      % look whether we got a fully dimensional polyhedron
                      H = [obj.Set(j).H; -obj.Set(i).H(k, :)];
                      if fast_isFullDim(H)
                          keep(k) = false;
                      end
                  end
              end
              keep_idx = find(keep);
              He = [He; obj.Set(i).H(keep_idx, :)];
          end
          
          if isempty(He)
              % the envelope is R^n
              E = Polyhedron.fullSpace(obj.Dim);
          else
              E = Polyhedron(He(:, 1:end-1), He(:, end));
          end
      end
      
      function result = compare(P1, P2, function_name)
          % Compares two PWA functions
          %
          %   result = P1.compare(P2, function_name)
          %
          % Output:
          %   0 if P1(x) == P2(x) for all x
          %   1 if P1(x) >= P2(x) for all x
          %   2 if P1(x) <= P2(x) for all x
          %   3 if P1(x) >= P2(x) for some x and P1(x) <= P2(x) for some
          %     other x
          
          global MPTOPTIONS
          
          narginchk(2, Inf);
          assert(isa(P1, 'PolyUnion'), 'The first input must be a PolyUnion object.');
          assert(isa(P2, 'PolyUnion'), 'The second input must be a PolyUnion object.');
          error(P1.rejectArray());
          error(P2.rejectArray());
          
          if nargin < 3 || isempty(function_name)
              % if no function is specified, take the first one
              if length(P1.listFunctions)>1
                  error('Please specify which function to use for comparison.');
              else
                  fnames = P1.listFunctions();
                  if isempty(fnames)
                      error('The object has no functions.');
                  end
                  function_name = fnames{1};
              end
          end
          
          if ~ischar(function_name)
			  error('The function name must be a string.');
		  elseif ~P1.hasFunction(function_name) || ~P2.hasFunction(function_name)
			  error('No such function "%s" in the object.', function_name);
		  elseif ~isa(P1.index_Set(1).Functions(function_name), 'AffFunction') || ...
                  ~isa(P2.index_Set(1).Functions(function_name), 'AffFunction')
			  error('Function "%s" must be affine.', function_name);
          elseif P1.index_Set(1).Functions(function_name).R~=1
              error('The function to compare must be scalar');
          elseif P1.index_Set(1).Functions(function_name).R ~= ...
              P2.index_Set(1).Functions(function_name).R
              error('Both functions must have the same range.');
          end

          if P1.Dim~=P2.Dim || P1.Domain~=P2.Domain
              error('The functions must be defined over the same domain.');
          end
          
          % since we only support PWA functions for now, it is enough to
          % compare the function values at the vertices
          vals = [];
          for i = 1:P1.Num
              V = P1.Set(i).V;
              for j = 1:size(V, 1)
                  x = V(j, :)';
                  fP1 = P1.feval(x, function_name, 'tiebreak', function_name);
                  fP2 = P2.feval(x, function_name, 'tiebreak', function_name);
                  vals = [vals; fP1(:) fP2(:)];
              end
          end          
          
          if norm(vals(:, 1)-vals(:, 2), Inf) < MPTOPTIONS.abs_tol
              % P1(x) = P2(x) for all x
              result = 0;
          elseif all(vals(:, 1) >= vals(:, 2))
              % P1(x) >= P2(x) for all x
              result = 1;
          elseif all(vals(:, 2) >= vals(:, 1))
              % P2(x) >= P1(x) for all x
              result = 2;
          else
              % P1(x)<=P2(x) and P2(y)<=P1(y)
              result = 3;
          end
      end
  end  
end

