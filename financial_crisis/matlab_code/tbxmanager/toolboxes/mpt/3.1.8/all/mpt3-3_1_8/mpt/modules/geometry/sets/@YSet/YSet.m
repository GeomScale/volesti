classdef YSet < ConvexSet
%
%  YSET: Representation of a convex set using YALMIP constraints. 
%  ===============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      S = YSet(vars, constraints)
%      S = YSet(vars, constraints, options)
%    
%  
%  DESCRIPTION
%  -----------
%     The class YSet represents convex sets described by YALMIP constraints.
%  Because YALMIP offers very broad specification of convex sets, the class YSet is
%  useful when applying methods of the ConvexSet class that are not available in
%  YALMIP. However, it is not intended with this class to replace basic
%  functionalities for YALMIP objects. For reference how to use YALMIP objects,
%  refer to YALMIP help. Only convex sets are accepted. Convexity is checked
%  internally by YALMIP.
%  
%  INPUT
%  -----
%     
%        
%          vars        Symbolic variables used in the           
%                      description of the constraints. The      
%                      dimension of the variables must much the 
%                      dimension used in the constraint set.    
%                      Vector and matrix variables are          
%                      accepted, multidimensional matrices are  
%                      not allowed.                             
%                      Class: sdpvar                            
%          constraints Constraint set given as lmi object. The  
%                      constraints must build a convex set,     
%                      otherwise the argument is not accepted.  
%                      The convexity is checked internally by   
%                      YALMIP.                                  
%                      Class: lmi                               
%          options     YALMIP options defined by sdpsettings.   
%                      You can specify the solver here,         
%                      verbosity, the tolerances, etc. By       
%                      default, these options are idependent of 
%                      MPT settings. YALMIP chooses the solver  
%                      based on its internal preferences and    
%                      depending on the type of the constraint  
%                      set. For more details, type help         
%                      sdpsettings.                             
%                      Class: struct                            
%                        
%  
%  
%  OUTPUT
%  ------
%     
%        
%          S YSet object representing a convex set.   
%            Class: YSet                              
%              
%  
%  
%  SEE ALSO
%  --------
%     ConvexSet,  Polyhedron
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
  properties(SetAccess = private)
    constraints; % YALMIP constraints
    vars;        % YALMIP variables
  end
  
  properties(SetAccess = public)
    opts = sdpsettings('verbose', 0,'allownonconvex',0);        % YALMIP default options
  end
  
  properties(SetAccess = private, Hidden = true)
    extr   = []; % YALMIP optimizer object for computing extreme points
    maxsep = []; % YALMIP optimizer object for projecting points
    alpha = sdpvar(1); % Used in shoot
    x0 = [];     % Used in extreme
  end
  
  %%
  methods

    %% Constructor
    function obj = YSet(vars, constraints, opts)
      %
      % obj = YSet(vars, constraints, opts)
      %
      % Constructor with default YALMIP options.<p>
      %
      % Note that the object can involve many YALMIP variables that are not listed in vars.
      % For example, one can operate directly on the projection of an object without first
      % projecting that object.
      %
      % @param vars The YALMIP variables involved in this object.
      % @param constraints The YALMIP constraints describing this object.
      %
      
	  if nargin==0
		  return
	  end
	  
      narginchk(2, 3);
      
      if ~isvector(vars)
          error('Variables must be provided as vectors only.');
      end
      if ~isa(vars, 'sdpvar'),
          error('YALMIP variables must be given as "sdpvar" object.');
      end

      if ~(isa(constraints, 'lmi') || isa(constraints, 'constraint'))
          error('YALMIP constraints must be given as "lmi" object.');
      end
      
      % assign arguments
      obj.constraints = constraints;
      obj.vars = vars;
	  % initialized the function storage
	  obj.Functions = containers.Map;
      
      if nargin > 2,
          if ~isa(opts,'struct')
              error('Options must be provided in a struct format, given by YALMIP "sdpsettings".');
          else
             % check the fields of sdpsettings
             fn = fieldnames(sdpsettings);
             for i=1:numel(fn)
                if ~isfield(opts,fn{i})
                    error('The field "%s" is missing in the options format.',fn{i});
                end
             end              
          end
          % for any supplied option we must always check for convexity
          obj.opts = sdpsettings(opts,'allownonconvex',0);
      end

      obj.Dim = numel(obj.vars);
      
      % check if the dimension matches with the dimension of the variables
      % in the constraint set
      
      % generate null point and assign
      v = zeros(size(obj.vars)); 
      s = yalmip('getsolution');
      assign(obj.vars,v);
      % check residuals
      r = checkset(obj.constraints);
      if any(isnan(r))
          error('Dimension mismatch between the provided variables and the variables in the constraint set.');
      end
      yalmip('setallsolution',s);
          

                  
      % Prep extr function
      % Create a model in YALMIPs low level format
      % All we change later is the cost vector
      [model,recoverdata,diagnostic,internalmodel] = export(obj.constraints,[],obj.opts,[],[],0);

      if ~isempty(diagnostic)
          % check convexity
          if diagnostic.problem==14
              error('Provided YALMIP constraints build non-convex set. Only convex set are allowed.');
          end
          % there may be another problem
          error('An error occured when exporting YALMIP model: "%s"', diagnostic.info);
      end

      if isempty(internalmodel)
        error('Could not create model for inner approximation.');
      end

      internalmodel.options.saveduals = 0;
      internalmodel.getsolvertime = 0;
      internalmodel.options.dimacs = 0;
      
      localindex = [];
      for i = 1:numel(obj.vars)
        localindex = [localindex find(ismember(recoverdata.used_variables,getvariables(obj.vars(i))))];
      end
      
      obj.extr.model = internalmodel;
      obj.extr.local = localindex;
      
      % Prep maxsep function
      obj.x0 = sdpvar(obj.Dim,1,'full');
      obj.maxsep = optimizer(obj.constraints, 0.5*(obj.x0-obj.vars(:))'*(obj.x0-obj.vars(:)),obj.opts,obj.x0,obj.vars);

      
	end
	        
  end
  
  methods(Access=protected)
	  
	  function new = copyElement(obj)
		% Copy method for YSet objects
		
		% Note: matlab.mixin.Copyable.copy() automatically properly
		% copies arrays and empty arrays, no need to do it here.

		% deep copy of Functions, shallow copy of Internal and Data
		new = copyElement@ConvexSet(obj);
		
		if numel(obj.vars)>0
			% create a copy of variables and constraints
			[new_vars, new_cons] = obj.copy_Y_constraints(obj.vars, obj.constraints);
		
			% use new variables and constraints
			new.vars = new_vars;
			new.constraints = new_cons;
		end
		
	end
        
  end
  
  methods (Static, Hidden)
	  
	  function [new_vars, new_con, new_con_s] = copy_Y_constraints(vars, constraints)
		  % Internal helper to create a copy of YALMIP's constraints
		  %
		  %   [new_vars, new_con] = mpt_copyConstraints(variables, constraints)
		  %
		  % Creates a copy of "constraints", which employ "variables". The output are
		  % new variables in "new_vars" and new constraints in "new_con"
		  %
		  % Examples:
		  %   x = sdpvar(2, 1);
		  %   C = [-1 <= x <= 1] + [ x'*x >= 2 ] + [ sin(cos(x(2))) <= x(1)^2 ];
		  %   [nx, nC, Cstr] = copy_Y_constraints(x, C)
		  %
		  %   x = sdpvar(2, 1); p = polytope(randn(5, 2), rand(5, 1));
		  %   C = [ ismember(x, p) ];
		  %   [nx, nC, Cstr] = copy_Y_constraints(x, C)
		  %
		  %   x = sdpvar(1, 1);
		  %   C = [-1 <= x <= 1] + [ x == 0 ];
		  %   [nx, nC, Cstr] = copy_Y_constraints(x, C)
		  
		  if size(vars, 2)~=1
			  error('Variable must be a column vector.');
		  end
		  
		  if is(vars, 'binary') || is(vars, 'integer')
			  error('Binary and integer variables not yet supported.');
		  end
		  
		  % create new variables
		  new_vars = sdpvar(length(vars), 1); % implicitly assumes a column vector
		  
		  % prepare string name of each component of "new_vars":
		  % 1) determine which YALMIP variables are part of the constraints
		  % 2) allocate enough "names"
		  % 3) set names of variables obtained in step 1
		  
		  % get total number of YALMIP sdpvar objects
		  s = yalmip('getinternalsdpvarstate');
		  n_total = numel(s.variabletype);
		  new_names = cell(1, n_total);
		  [new_names{:}] = deal('');
		  
		  % set names of variables of interest
		  vars_idx = getvariables(vars);
		  local_idx = 0;
		  for i = vars_idx
			  local_idx = local_idx + 1;
			  % implicitly assumes the variables are contiguous
			  % TODO: deal with non-contiguous variables
			  new_names{i} = sprintf('new_vars(%d)', local_idx);
		  end
		  
		  % separate equalities from inequalities
		  idx_eq = find(is(constraints, 'equality'));
		  idx_ineq = setdiff(1:length(constraints), idx_eq);
		  new_con = [];
		  new_con_s = cell(length(sdpvar(constraints)), 1);
		  old_cons = { constraints(idx_eq), constraints(idx_ineq) };
		  signs = { '==', '>=' };
		  idx = 1;
		  for i = 1:numel(old_cons)
			  % create new constraints
			  if length(old_cons{i})==0, continue, end
			  con_s = sdisplay2(sdpvar(old_cons{i}), new_names);
			  if ~iscell(con_s)
				  % just for future. right now the output of sdisplay2() is
				  % already a cell
				  con_s = { con_s };
			  end
			  for j = 1:numel(con_s)
				  new_con_s{idx} = sprintf('%s %s 0', con_s{j}, signs{i});
				  new_con = new_con + eval(new_con_s{idx});
				  idx = idx + 1;
			  end
		  end
		  
	  end
 
  end
  
end

