classdef EMPCController < AbstractController
%
%  EMPCCONTROLLER: Explicit MPC controller 
%  ========================================
%  
%  
%  SYNTAX
%  ------
%     
%      ctrl = EMPCController(model, horizon)
%      ctrl = EMPCController(MPCController)
%    
%  
%  DESCRIPTION
%  -----------
%     Constructs the explicit form of an MPC controller. The particular type of the
%  optimization problem to be solved parametrically is determined by the type of
%  the prediction model and by its parameters. For a more detailed information, see
%  " help MPCController".
%    Instances of the EMPCController class expose following properties: 
%    
%     - model: the prediction model used in the MPC setup; 
%     - N: the prediction horizon 
%     - optimizer: the explicit optimizer as an instance of the PolyUnion class; 
%     - partition: the polyhedral partition of the explicit feedback law as an
%     instance of the Polyhedron class; 
%     - feedback: the explicit representation of the feedback law as an instance of
%     the PolyUnion class; 
%     - cost: the explicit representation of the optimal cost function as an
%     instance of the PolyUnion class. 
%    The optimizer property is available for read/write access. This allows, for
%  instance, to remove overlaps from multiple overlapping partitions by
%  ctrl.optimizer = ctrl.optimizer.merge().
%  
%  INPUT
%  -----
%     
%        
%          model Any MPT3 system ( LTISystem, PWASystem,  
%                MLDSystem)                               
%                Class: AbstractSystem                    
%                  
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ctrl Explicit MPC controller                  
%                 
%  
%  
%  SEE ALSO
%  --------
%     MPCController
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
 
 
	%   ctrl = EMPCController(model, N)
	%   ctrl = EMPCController(MPCController)
    properties(SetAccess=protected, Transient=true)
		partition
		feedback
		cost
		nr
	end
    
    methods

		function out = isExplicit(obj)
            % True for explicit controllers
            out = true;
		end
		
		function out = getName(obj)
			out = 'Explicit MPC controller';
		end

		function obj = EMPCController(varargin)
            % Constructor:
			%
			%   ctrl = EMPCController(model, N)

			obj.addlistener('optimizer', 'PostSet', @obj.optPostSetEvent);

			if nargin==0
                return
			end

			obj.importUserData(varargin{:});
			
			if ~isobject(obj.optimizer)
				obj.construct();
			end
		end
		
		function obj = toInvariant(in)
			% Computes the invariant subset of the controller
			
			if nargout==0
				% replace original controller
				obj = in;
			else
				% create a new controller
				obj = in.copy();
			end
			
			% create the closed-loop model
			CL = ClosedLoop(obj, obj.model).toSystem();
			
			% compute invariant subset
			[I, dynamics] = CL.invariantSet();
			
			% create new explicit solution
			for j = 1:length(obj.optimizer.Set(1).Func)
				for i = 1:length(I)
					I(i).copyFunctionsFrom(obj.optimizer.Set(dynamics(i)));
				end
			end
			
			% TODO: propagate properties? (Overlaps, Bounded, FullDim)
			obj.optimizer = PolyUnion('Set', I);
			
		end
		
        function obj = construct(obj)
            % Constructs the explicit solution

			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			% make sure the prediction horizon was provided
			error(obj.assert_controllerparams_defined);

			if isempty(obj.yalmipData)
				Y = obj.toYALMIP();
			else
				Y = obj.yalmipData;
			end
			% TODO: remove the debug flag for beta
			options = sdpsettings('debug', 1);

			speedhack = false;

			if speedhack
				% bypass Opt() since it's so slow now (Aug 2012)
				sol = solvemp(Y.constraints, Y.objective, options, ...
					Y.internal.parameters, Y.variables.u(:));
				if isempty(sol) || isempty(sol{1})
					error('Problem is infeasible.');
				end
				% convert to MPT3
				res.xopt = mpt_mpsol2pu(sol);
			else
				% standard way via Opt

				% TODO: Opt() should allow MILP/MIQP problems to be constructed
				% even if MPLP/MPQP parametric solvers are chosen
				problem = Opt(Y.constraints, Y.objective, ...
					Y.internal.parameters, Y.variables.u(:));
				
				if isequal(problem.problem_type,'MILP') || ...
						isequal(problem.problem_type, 'MIQP')
					% use YALMIP for pMILP/pMIQP problems
					sol = solvemp(Y.constraints, Y.objective, options, ...
						Y.internal.parameters, Y.variables.u(:));
					if isempty(sol) || isempty(sol{1})
						error('Problem is infeasible.');
					end
					% convert to MPT3
					res.xopt = mpt_mpsol2pu(sol);
				else
					res = problem.solve;
					if res.exitflag == MPTOPTIONS.INFEASIBLE
						error('Problem is infeasible.');
					end
				end
			end
            obj.optimizer = res.xopt;
		end
        
		% 		function h = plot(obj)
		% 			% Plots partition of the controller
		%
		% 			h = obj.optimizer.plot();
		% 			if nargout==0
		% 				clear h
		% 			end
		% 		end
		
		function [u, feasible, openloop] = evaluate(obj, xinit, varargin)
            % Evaluates the explicit solution for a given point
			%
			% u = controller.evaluate(x0) evaluates the explicit MPC
			% solution and returns the optimal control input associated to
			% point "x0". If "x0" is outside of the controller's domain,
			% "u" will be NaN.
			%
			% [u, feasible] = controller.evaluate(x0) also returns a
			% feasibility boolean flag.
			%
			% [u, feasible, openloop] = controller.evaluate(x0) also
			% returns the full open-loop optimizer in "openloop.U" and the
			% optimal cost in "openloop.cost". Moreover, "openloop.X" and
			% "openloop.Y" will be set to NaN. This is due to the fact that
			% storing open-loop predictions of states and outputs in the
			% explicit solution would considerably increase its size.
			%
			% u = controller.evaluate(x0, 'x.reference', ref, 'u.prev', u0)
			% includes "ref" and "u0" into the vector of initial
			% conditions.

			% TODO: point location should be a method of PolyUnion

			% make sure the prediction horizon was provided
			error(obj.assert_controllerparams_defined);

			if isempty(obj.nu)
				error('Set ctrl.nu first.');
			end
			
            % evaluate the explicit optimizer
			error(validate_vector(xinit, obj.nx, 'initial state'));
			
			% assemble the vector of initial conditions. Include any
			% variables that were declared by filters as those which need
			% proper initialization.
			xinit = obj.parse_xinit(xinit, varargin{:});
			
			% evaluate the primal optimizer, break ties based on the cost
			% function. guarantees that the output is a single region where
			% the cost is minimal.
			[U, feasible, idx, J] = obj.optimizer.feval(xinit, ...
				'primal', 'tiebreak', 'obj');
			if ~feasible
				J = Inf;
				% index of the optimizer and index of the region from which the
				% control action was extracted
				opt_partition = [];
				opt_region = [];

			elseif numel(obj.optimizer)==1
				opt_partition = 1;
				opt_region = idx;

			else
				% multiple optimizers
				opt_partition = idx(1);
				opt_region = idx(2);
			end
				
			if isempty(J)
				% no tie-breaking was performed, compute cost manually
				%
				% Note: from a long-term sustainibility point of view
				% we should use
				%   J = obj.optimizer.Set(idx).feval(xinit, 'obj');
				% here. but ConvexSet/feval() adds so much unnecessary
				% overhead that we better evaluate the function
				% directly
				J = obj.optimizer(opt_partition).Set(opt_region).Functions('obj').feval(xinit);
			end				
			
			if numel(U)~=obj.nu*obj.N
				% sanity checks for EMPCControllers imported from
				% polyunions
				error('Number of optimizers is inconsistent with "N" and/or "nu".');
			end
			
			u = U(1:obj.nu);
			if nargout==3
				openloop.cost = J;
				openloop.U = reshape(U, [obj.nu obj.N]);
				openloop.X = NaN(obj.nx, obj.N+1);
				openloop.Y = NaN(obj.model.ny, obj.N);
				openloop.partition = opt_partition;
				openloop.region = opt_region;
			end
		end
		
		function new = simplify(obj, method)
			% Simplifies a given explicit controller
			%
			%   simple = controller.simplify(method)
			%
			% Supported methods:
			%    'clipping': Clipping approach
			%                (see ClippingController)
			%     'fitting': PWA fitting controller
			%                (see FittingController)
			%      'greedy': Greedy region merging
			%         'orm': Optimal Region Merging
			%  'separation': Saturation-based separation
			%                (see SeparationController)
			
			narginchk(1, Inf);
			if nargin < 2
				method = 'greedy';
			end
			
			% trim the feedback law to just the receding horion part.
			% operate on a copy as not to destroy the original object
			old = obj.copy();
			old.optimizer.trimFunction('primal', obj.nu);
			old.N = 1; % to get correct size of the open-loop optimizer
			% TODO: implement a better way
			
			switch lower(method)
				case {'greedy', 'orm'},
					optimal = isequal(lower(method), 'orm');
					if nargout==0
						new = obj;
					else
						new = old;
					end
					new.optimizer = old.optimizer.merge('primal', 'optimal', optimal);
					
					% TODO: figure out why the event is not triggered
					% automatically
					new.optPostSetEvent();
					
				case 'fitting',
					new = FittingController(old);
					
				case 'separation',
					new = SeparationController(old);
					
				case 'clipping',
					new = ClippingController(old);
					
				otherwise,
					help EMPCController/simplify
					error('Unrecognized simplification method "%s".', method);
			end
		end
		
		function data = clicksim(obj, varargin)
			% Select initia condition for closed-loop simulation by mouse
			%
			%   controller.clicksim(['option', value, ...)
			%
			% Select initial points by left-click. Abort by right-click.
			%
			% options:
			%  'N_sim': length of the closed-loop simulation (default: 100)
			%  'x0': initial point, if provided, the method exits
			%        immediately (default: [])
			%  'model': model for the closed-loop simulation
			%           (default is controller.model)
			%  'alpha': transparency of the partition (default: 1)
			%  'color': color of the closed-loop profiles (default: 'k')
			%  'linewidth': width of the line (default: 2)
			%  'marker': markers indicating points (default: '.')
			%  'markersize': size of the markers (default: 20)
			
			if obj.nx~=2
				error('Only 2D partitions can be plotted.');
            end
            
            % R2016b specific: inputParser does not accept 'x.something'
            inputs = varargin;
            for i = 1:2:length(inputs)
                % 'x.something' -> 'x_something'
                % but we only use "inputs" for option parsing and we still
                % feed the original "varargin" to ClosedLoop/simulate()
                inputs{i} = strrep(inputs{i}, '.', '_');
            end
            
			ip = inputParser;
			% important: use keepUnmatched=true to propage free references
			% to ctrl.simulate()
			ip.KeepUnmatched = true;
			ip.addParamValue('x0', []);
			ip.addParamValue('N_sim', 50, @isnumeric);
			ip.addParamValue('alpha', 1, @isnumeric);
			ip.addParamValue('linewidth', 2, @isnumeric);
			ip.addParamValue('color', 'k');
			ip.addParamValue('marker', '.');
			ip.addParamValue('markersize', 20);
			ip.addParamValue('model', obj.model);
            ip.addParamValue('openloop', false);
			ip.parse(inputs{:});
			options = ip.Results;

			% closed-loop system
			loop = ClosedLoop(obj, options.model);
			
			% if we have free references or other variables which extend
			% the vector of initial conditions, plot the section
			if isstruct(obj.xinitFormat) && obj.xinitFormat.n_xinit>obj.nx
				P = obj.optimizer.Set.slice(obj.nx+1:obj.optimizer.Dim);
			else
				P = obj.partition;
			end
			P.plot('alpha', options.alpha);
			axis tight
			hold on
			button = 1;
			while button~=3
				if isempty(options.x0)
					[x, y, button] = ginput(1);
					x0 = [x; y];
				else
					x0 = options.x0;
					button = 3;
                end
                if options.openloop
                    [~, ~, ol] = obj.evaluate(x0);
                    model = options.model.copy();
                    model.initialize(x0);
                    data = model.simulate(ol.U);
                else
                    data = loop.simulate(x0, options.N_sim, varargin{:});
                end
				plot(data.X(1, :), data.X(2, :), options.color, ...
					'linewidth', options.linewidth, ...
					'marker', options.marker, ...
					'markersize', options.markersize);
				title(sprintf('Closed-loop simulation: x0 = %s', mat2str(x0)));
			end
			hold off
			if nargout==0
				clear data
			end
			
		end
		
		function out = binaryTree(obj)
			% Creates a binary tree representation of the optimizer
			%
			% ctrl.binaryTree() converts the controller's optimizer into a
			% BinTreePolyUnion object. Subsequent evaluation via
			% ctrl.evaluate() will then automatically use
			% BinTreePolyUnion/contains() to perform point location.
			%
			% new = ctrl.binaryTree() will leave "ctrl" intact, and instead
			% creates a copy "new", and converts that one into binary tree.
			
			if nargout>0
				out = obj.copy();
			else
				out = obj;
			end
			newopt = [];
			for i = 1:numel(obj.optimizer)
				newopt = [newopt, BinTreePolyUnion(obj.optimizer(i))];
			end
			out.optimizer = newopt;
		end
		
	end
	
	methods(Static, Hidden)
		% private APIs, use at your own risk
		
		function optimizer = addMissingFunctions(optimizer, varargin)
			% Adds missing functions to the optimizer
			
			% TODO: check dimensions
			
			% Any PolyUnion representing an explicit solution has to
			% contain following functions: 'z', 'w', 'primal', 'dual',
			% 'obj'. Any missing function is replaced by a zero function.
			%
			% Moreover, we can add other functions, specified in the
			% function call, e.g.:
			%    opt = obj.addMissingFunctions(obj, opt, 'myFun',
			%    AffFunction)
			% which adds the function 'myfun' to all polyhedra of the
			% PolyUnion 'opt'.

			required_names = {'z', 'w', 'primal', 'dual', 'obj'};
			required_args = cell(1, length(required_names));
			for i = 1:2:length(varargin)/2
				required_names{end+1} = varargin{i};
				required_args{end+1} = varargin{i+1};
			end
			
			% TODO: keep the list of functions in sync with the Opt class
			for i = 1:length(optimizer)
				funs = optimizer(i).listFunctions;
				missing = find(ismember(required_names, funs)==0);
				if ~isempty(missing)
					for k = 1:length(missing)
						for j = 1:length(optimizer(i).Set)
							if isempty(required_args{missing(k)})
								addf = AffFunction(zeros(1, optimizer(i).Set(j).Dim));
							else
								addf = required_args{missing(k)};
							end
							optimizer(i).Set(j).addFunction(addf, ...
								required_names{missing(k)});
						end
					end

					% update list of functions associated to the union
					optimizer(i).setInternal('FuncName', required_names);
				end
			end
		end
		
		function [nx, nu] = checkOptimizer(optimizer)
			% checks sanity of an explicit optimizer
			
			if isa(optimizer, 'double') && isempty(optimizer)
				% setting the optimizer to [] should work
				return
			end
			if ~isa(optimizer, 'PolyUnion')
				error('Optimizer must be an instance of the @PolyUnion class.');
			end
			
			nx = zeros(1, numel(optimizer));
			nu = zeros(1, numel(optimizer));
			% check that we have at least the 'primal' and 'obj'
			% functions and determine number of paramters (nx) and number
			% of optimizer's outputs (nu)
			for i = 1:numel(optimizer)
				if optimizer(i).Num<1
					error('Optimizer %d must contain at least one region.', i);
				elseif ~any(cellfun(@(x) isequal(x, 'primal'), optimizer(i).listFunctions))
					error('Optimizer %d must contain the "primal" function.', i);
				elseif ~any(cellfun(@(x) isequal(x, 'obj'), optimizer(i).listFunctions))
					error('Optimizer %d must contain the "obj" function.', i);
				end
				nx(i) = optimizer(i).Dim;
				primal = optimizer(i).Set(1).getFunction('primal');
				nu(i) = length(primal.g);
			end
			% check that all polyunions are of the same dimension
			if any(diff(nx)~=0)
				error('All optimizers must have equal dimensions.');
			end
			if any(diff(nu)~=0)
				error('All optimizers must have the same number of optimization variables.');
			end
			nx = nx(end);
			nu = nu(end);
		end

	end
	
	methods(Access=protected)
		
		function obj = importUserData(obj, varargin)
			% Imports model / prediction horizon

			obj = obj.importUserData@AbstractController(varargin{:});

			if nargin==2 
				% possible import from an another controller
				otherCtrl = varargin{1};
				if isa(otherCtrl, 'AbstractController') && ...
						otherCtrl.isExplicit() && isobject(otherCtrl.optimizer)
					% import optimizer
					obj.model = otherCtrl.model;
					obj.N = otherCtrl.N;
					obj.optimizer = otherCtrl.optimizer;
					%obj.optimizer=obj.addMissingFunctions(otherCtrl.optimizer);
					return
					
				elseif isa(otherCtrl, 'PolyUnion')
					% import optimizer, use dummy system
					obj.model = LTISystem;
					optimizer = otherCtrl;
					obj.nx = obj.checkOptimizer(optimizer);
					obj.optimizer = optimizer.copy;
					%obj.optimizer=obj.addMissingFunctions(otherCtrl);
					return
				end
			end

		end

	end

	methods(Access=protected, Hidden)
	
		function optPostSetEvent(obj, source, event)
			% triggered after the optimizer was changed
			%
			% store the partition, feedback and cost objects for fast
			% access

			
			obj.markAsModified();

			if isempty(obj.optimizer)
				obj.feedback = [];
				obj.cost = [];
				obj.feedback = [];
				return
			end
			
			obj.feedback = obj.optimizer.getFunction('primal');
			obj.cost = obj.optimizer.getFunction('obj');
			if numel(obj.optimizer)==1
				% single optimizer, copy it and remove functions
				out = obj.optimizer.copy();
				out.removeAllFunctions();
			else
				% multiple optimizers, concatenate regions together
				P = [];
				for i = 1:length(obj.optimizer)
					P = [P; Polyhedron(obj.optimizer(i).Set)];
				end
				out = PolyUnion('Set', P.removeAllFunctions, ...
					'Domain', cat(1, obj.optimizer.Domain));
			end
			obj.partition = out;
			
			% number of regions
			if isobject(obj.optimizer)
				n = cell(1, numel(obj.optimizer));
				[n{:}] = obj.optimizer.Num;
				obj.nr = sum([n{:}]);
			else
				obj.nr = 0;
			end
		end

	end
	
	methods(Static)
		
		function new = loadobj(obj)
			% load method
			
			% post-set events must be triggered manually
			new = obj;
			new.optPostSetEvent();
		end

	end
end
