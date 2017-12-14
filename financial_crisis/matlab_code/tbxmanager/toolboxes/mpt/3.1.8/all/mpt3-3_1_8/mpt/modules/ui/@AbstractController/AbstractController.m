classdef AbstractController < FilterBehavior & ComponentBehavior & IterableBehavior
    % Class representing an abstract controller
    %
    % Each derived class has to define at least following methods:
    %   - constructor
    %   - evaluate() which returns the optimizer
    %   - isExplicit() which returns true for controllers which can be
    %                  converted to C code

	properties(AbortSet, SetObservable=true, GetObservable=true)
		N % Prediction horizon
	end
	properties(SetObservable=true, GetObservable=true)
		model % Prediction model
		optimizer % Implicit or explicit optimizer
	end
	properties
		nu % Dimension of the optimizer
		nx % Dimension of the paramater
	end
	properties(SetAccess=protected, Hidden)
		yalmipData % YALMIP representation of the MPC problem
		xinitFormat % format of the initial condition
	end

	methods(Abstract)
		% List of methods which any controller has to implement
		
		% compute control action for a given value of the initial condition
		[u, feasible, openloop] = evaluate(obj, xinit)
		
		% returns true if controller is explicit
		out = isExplicit(obj)
		
		% return human-readable name of the control scheme
		out = getName(obj)
		
		% constructs the implicit or explicit optimizer
		obj = construct(obj)
		
		% validate the optimizer
		checkOptimizer(optimizer)
	end

	methods

		function obj = AbstractController(varargin)
			% Constructor
        end
        
        function out = mtimes(A, B)
            % A*B creates a closed-loop system
            
            if isa(A, 'AbstractSystem') && isa(B, 'AbstractController')
                model = A;
                ctrl = B;
            elseif isa(A, 'AbstractSystem') && isa(B, 'double')
                model = A;
                ctrl = SFController(model, B);
            elseif isa(B, 'AbstractSystem') && isa(A, 'AbstractController')
                model = B;
                ctrl = A;
            elseif isa(B, 'AbstractSystem') && isa(A, 'double')
                model = B;
                ctrl = SFController(model, A);
            else
                error('Unsupported inputs.');
            end
            out = ClosedLoop(ctrl, model);
        end

        function obj = set.optimizer(obj, value)
			% validate and set the optimizer
		
			obj.checkOptimizer(value);
			obj.optimizer = value;
			obj.optPostSetEvent();
		end
		
		function display(obj)
			% Display
			
			if numel(obj)>1
				fprintf('Array of %d %ss\n', numel(obj), obj.getName());
			elseif isempty(obj.model)
				fprintf('Empty %s\n', obj.getName());
			else
				fprintf('%s', obj.getName());
				if ~isempty(obj.N)
					fprintf(' (horizon: %d', obj.N);
					if obj.isExplicit()
						if obj.nr>0
							fprintf(', regions: %d', obj.nr);
							if numel(obj.optimizer)>1
								fprintf(', partitions: %d', numel(obj.optimizer));
							end
						else
							fprintf(', explicit solution not yet available');
						end
					end
					fprintf(')\n');
				else
					fprintf(' (no prediction horizon defined)\n');
				end
			end

		end
		
		function obj = set.model(obj, model)
			% Sets the prediction model
			
			obj.model = model.copy();
			obj.model.prepareForControl();
			obj.nu = model.nu;
			obj.nx = model.nx;
		end
		
		function obj = set.N(obj, N)
			% Validates value of the prediction horizon
			
			error(obj.validatePredictionHorizon(N));
			obj.N = N;
		end			


		function Y = toYALMIP(obj)
			% Converts the optimization setup into YALMIP's objects
			%
			% NOTE: derived classes which implement control algorithms that
			% do not have a straightforward YALMIP implementation should
			% re-implement this method to throw an error
			
			% make sure we have the prediction horizon available
			error(obj.assert_controllerparams_defined);
			
			Y = struct('constraints', [], 'objective', [], ...
				'variables', [], 'internal', []);
			obj.model.instantiate(obj.N);
			Y.constraints = obj.model.constraints();
			Y.objective = obj.model.objective();
			
			% TODO: automatically include other variables using
			% getComponents
			Y.variables = struct('x', obj.model.getVariables('x'), ...
				'u', obj.model.getVariables('u'), ...
				'y', obj.model.getVariables('y'));
			
			if isa(obj.model, 'PWASystem')
				Y.variables.d = obj.model.getVariables('d');
				
			elseif isa(obj.model, 'MLDSystem')
				Y.variables.d = obj.model.getVariables('d');
				Y.variables.z = obj.model.getVariables('z');
			end
			
			% list of variables which we ask for from the optimizer in the
			% following order:
			%  * cost (scalar)
			%  * inputs (nu*N)
			%  * states (nx*N)
			%  * outputs (ny*N)
			%  * all other model variables
			f = fieldnames(Y.variables);
			% the cost must always be the first variable
			vars = Y.objective;
			main_variables = {'u', 'x', 'y'};
			for i = 1:length(main_variables)
				v = Y.variables.(main_variables{i});
				vars = [vars; v(:)];
			end
			
			% now add all other variables (e.g. 'd' for PWA, 'd','z' for
			% MLD), dive recursively into automatically introduced
			% variables
			%
			% TODO: no need to include all variables since we only need
			% "u", "x", "y" in evaluate()
			for i = 1:length(f)
				if ~ismember(f{i}, main_variables)
					v = Y.variables.(f{i});
					vars = struct2vars(v(:), vars);
				end
			end
			
			% include additional variables introduced by filters
			add = containers.Map;
			keys = fields(Y.variables);
			for i = 1:length(keys)
				add(keys{i}) = obj.model.(keys{i}).applyFilters('getVariables', 'map');
			end
			% always add 'model' filters
			add('model') = obj.model.applyFilters('getVariables', 'map');
			% any variables introduced by filters?
			new_variables = false;
			keys = add.keys;
			for i = 1:length(keys)
				if ~isempty(add(keys{i}))
					new_variables = true;
					break
				end
			end
			
			% initial conditions for the optimization
			init_vars = struct('component', 'x', ...
				'name', 'x.init', ...
				'var', obj.model.x.var(:, 1), ...
				'dims', size(obj.model.x.var(:, 1)));
			
			if new_variables
				% include variables introduced by filters
				[filter_vars, init_vars] = map2struct(add, init_vars);
				Y.variables.filters = filter_vars;
			end
			
			% store format of the initial condition
			xinit_variables = [];
			obj.xinitFormat.names = {};
			obj.xinitFormat.components = {};
			obj.xinitFormat.dims = {};
			for i = 1:length(init_vars)
				xinit_variables = [xinit_variables; init_vars(i).var(:)];
				obj.xinitFormat.names{end+1} = init_vars(i).name;
				obj.xinitFormat.components{end+1} = init_vars(i).component;
				obj.xinitFormat.dims{end+1} = init_vars(i).dims;
			end
			% number of required initial conditions
			obj.xinitFormat.n_xinit = length(xinit_variables);

			% sdpvars of initial conditions
			Y.internal.parameters = xinit_variables;
			% sdpvars of all requested variables
			Y.internal.requested = vars;
			% format of the initial condition
			Y.internal.xinitFormat = obj.xinitFormat;
            % default sdpsettings
            Y.internal.sdpsettings = sdpsettings('verbose', 0);

		end
		
		function obj = fromYALMIP(obj, Y)
			% Import a YALMIP representation of constraints and objective
			%
			% NOTE: derived classes which implement control algorithms that
			% do not have a straightforward YALMIP implementation should
			% re-implement this method to throw an error

			if ~isstruct(Y) || ~isfield(Y, 'constraints') || ...
					~isfield(Y, 'objective') || ...
					~isfield(Y, 'variables') || ...
					~isfield(Y, 'internal')
				error('Input must be a structure generated by toYALMIP.');
			end
			% update the objective
			Y.internal.requested(1) = Y.objective;
			
			obj.xinitFormat = Y.internal.xinitFormat;
			obj.yalmipData = Y;
			obj.construct();
		end
		
		function new = saveobj(obj)
			% save method

			% warn about properties which cannot be saved
			new = copy(obj);
			if ~isempty(obj.yalmipData)
				fprintf('Custom YALMIP setups cannot be saved.\n');
				new.yalmipData = [];
			end
			
			new.saveAllFilters();
			new.saveAllComponents();
			
			% NOTE: obj.model will be properly saved by
			% AbstractSystem/saveobj

		end
		
		function out = simulate(obj, x0, N_sim, varargin)
			% Simulate the closed-loop system using the prediction model
			%
			%   data = controller.simulate(x0, N_sim)
			%
			% Inputs:
			%   controller: a controller object
			%   x0: initial point for the simulation
			%   N_sim: number of simulation steps
			%
			% Outputs:
			%   data: structure containing closed-loop profiles of states,
			%         inputs, outputs, and the cost function
			%
			% See ClosedLoop/simulate for more information.
			
			out = ClosedLoop(obj, obj.model).simulate(x0, N_sim, varargin{:});
        end
        
        function [avg_p, max_p, P1, P2, X0] = comparePerformance(this, other, varargin)
            % Compares performance of two controllers
            %
            %   [avg_diff, max_diff] = ctrl.comparePerformance(other)
            %
            % Compares the performance of two controllers by simulating
            % their closed-loop response. The simulations are performed for
            % a grid of initial conditions. Each closed-loop response is
            % then evaluated as \sum x_i'*Qx*x_i+u_i'*Qu*u_i+y_i'*Qy*y_i.
            %
            % The function returns five outputs:
            %   avg_diff: average relative difference of performance (in %)
            %   max_diff: worst-case relative difference (in %)
            %         P1: performance of the baseline controller (each entry
            %             indicates the closed-loop performance for the
            %             corresponding initial condition)
            %         P2: performance of the other controller (each entry
            %             indicates the closed-loop performance for the
            %             corresponding initial condition)
            %         X0: list of initial conditions
            %
            % Various options can be specified by 
            %   ctrl.comparePerformance(other, 'option1', value1, ...)
            %
            % List of options:
            %       'X0': array of initial conditions to investigate
            %             (if provided, the "domain" and "grid" options are
            %             ignored)
            %   'domain': the range of initial conditions to investigate as
            %             a Polyhedron object (if one of the controllers is
            %             an explicit one, the domain is taken from
            %             ctrl.optimizer.Domain)
            %     'grid': density of the grid of initial conditions 
            %             (default is 20)
            %     'Nsim': number of simulation steps
            %             (default is 20)
            %       'Qx': penalty matrix to evaluate performance
            %             (default is eye(nx))
            %       'Qu': penalty matrix to evaluate performance
            %             (default is eye(nu))
            %       'Qy': penalty matrix to evaluate performance
            %             (default is zeros(ny))
            
            assert(isa(other, 'AbstractController'), 'The other controller must be an MPT controller.');
            
            ip = inputParser;
            ip.addParamValue('domain', [], @validate_polyhedron);
            ip.addParamValue('grid', 10, @isscalar);
            ip.addParamValue('Nsim', 20, @isscalar);
            ip.addParamValue('Qx', eye(this.model.nx), @validate_realmatrix);
            ip.addParamValue('Qu', eye(this.model.nu), @validate_realmatrix);
            ip.addParamValue('Qy', zeros(this.model.ny), @validate_realmatrix);
            ip.addParamValue('X0', [], @validate_realmatrix);
            ip.parse(varargin{:});
            options = ip.Results;
            if isempty(options.X0)
                if isempty(options.domain)
                    if isa(this, 'EMPCController')
                        options.domain = this.optimizer.Domain;
                    elseif isa(other, 'EMPCController')
                        options.domain = other.optimizer.Domain;
                    else
                        error('The "domain" parameter must be specified.');
                    end
                end
                options.X0 = options.domain.grid(options.grid);
            end
            
            % performance function
            perf = @(d) sum(diag(d.X'*options.Qx*d.X)) + ...
                sum(diag(d.U'*options.Qu*d.U)) + ...
                sum(diag(d.Y'*options.Qy*d.Y));
            
            % evaluate performance for each starting point
            nx0 = size(options.X0, 1);
            P1 = zeros(nx0, 1);
            P2 = zeros(nx0, 1);
            for i = 1:nx0
                if i==1
                    fprintf('Comparing performance  ...');
                else
                    fprintf(1, ' \b\b\b\b\b%2d %%', fix(i/nx0*100));
                end
                if i==nx0
                    fprintf('\n');
                end
                x0 = options.X0(i, :)';
                d1 = this.simulate(x0, options.Nsim);
                d2 = other.simulate(x0, options.Nsim);
                P1(i) = perf(d1);
                P2(i) = perf(d2);
            end
            
            % round to 10 decimal places to avoid numerical noise
            P1 = round(P1*1e10)/1e10;
            P2 = round(P2*1e10)/1e10;
            
            d = abs(P1-P2)./P1;
            % avoid division by zero
            d(isnan(d)) = 0;
            avg_p = mean(d)*100;
            max_p = max(d)*100;
            
            X0 = options.X0;
        end


	end

	methods(Static, Hidden)

		function U = PolyhedronToPolyUnionHelper(P)
			% Helper function to convert a polyhedron into a
			% control-compatible polyunion.

			% After changeset 3f21e875f2ea, all polyunions imported to
			% EMPCController need to define at least the 'primal' and 'obj'
			% functions. This helper assigns these functions to an input
			% polyhedron 'P' (which can also be an array of polyhedra) and
			% returns a polyunion 'U'.

			% make a copy before adding functions
			P = Polyhedron(P);

			% add the functions
			for i = 1:numel(P)
				P(i).addFunction(AffFunction(zeros(1, P(i).Dim), 0), 'primal');
				P(i).addFunction(AffFunction(zeros(1, P(i).Dim), 0), 'obj');
			end

			U = PolyUnion(P);
		end

		function msg = validatePredictionHorizon(N, varargin)
			% Validates value of the prediction horizon
			
			ip = inputParser;
			ip.addParamValue('allowInf', true, @islogical);
			ip.parse(varargin{:});
			options = ip.Results;
			
			if isinf(N) && ~options.allowInf
				msg = 'The prediction horizon must be finite.';
			elseif N<1
				msg = 'The prediction horizon must be positive.';
			elseif round(N)~=N
				msg = 'The prediction horizon must be an integer.';
			else
				msg = '';
			end
		end
		
	end

	methods(Access=protected)

		function new = copyElement(obj)
			% Copy constructor

			new = copyElement@FilterBehavior(obj);
			
			% either copy the optimizer here, or move copying to
			% set.optimizer()
			if ismethod(obj.optimizer, 'copy')
				new.optimizer = obj.optimizer.copy();
			end

			% make a deep copy of obj.model
			%
			% Note: the set.model() method already performs a deep copy
			% when the model is changed, so the following line could be
			% changed to
			%    new.model = obj.model;
			% However, in terms of long-term sustainibility of the code we
			% better do the copy here as well, just in case we later remove
			% it from set.model().
			%
			% Another option would be to add the 'model' property
			% dynamically using the "addComponent()" method. In such a case
			% we could remove the whole AbstractController/copyElement()
			% method, since deep copy of all components is performed via
			% MPTUIHandle/copyElement()
			
			new.model = copy(obj.model);
		end

		function obj = importUserData(obj, model, predictionHorizon)
			% Imports model / prediction horizon
			
			if nargin==1
				return
			end
			
			if nargin==2 && isa(model, 'AbstractController')
				% import from an another controller
				otherCtrl = model;
				obj.model = otherCtrl.model;
				obj.N = otherCtrl.N;
				obj.xinitFormat = otherCtrl.xinitFormat;
				
			elseif isa(model, 'AbstractSystem')
				% import model and possibly the prediction horizon
				%obj.N = 1;
				obj.model = model;
				if obj.model.nu<0
					error('Cannot control autonomous systems.');
				end
				if nargin==3
					if isinf(predictionHorizon)
						error('The prediction horizon must be finite.');
					end
					obj.N = predictionHorizon;
				end
			end
		end
		
		function msg = assert_controllerparams_defined(obj)
			% Raises an error if the prediction horizon and/or other
			% parameters (number of states, number of inputs) are not
			% specified
			
			if isempty(obj.N) || obj.N==0
				msg = 'The prediction horizon must be specified.';
			elseif isempty(obj.nu)
				msg = 'Number of inputs must be specified.';
			elseif isempty(obj.nx)
				msg = 'Number of states must be specified.';
			else
				msg = '';
			end
			if nargout==0
				% no outputs = directly trigger the error
				% one output = just return the message
				error(msg);
			end
		end
	end

	methods(Access=protected, Hidden)
	
		function optPostSetEvent(obj, source, event)
			% triggered after the optimizer was changed
			
			obj.markAsModified();
		end
	end

	methods(Hidden)
		% internal APIs
		
		function xinit = parse_xinit(obj, x0, varargin)
			% Constructs initial conditions for the optimizer
			%
			%    xinit = ctrl.parse_xinit(x0, 'u.prev', u0, 'x.reference', r)
			%
			% Returns xinit = { x0, r, u0 }, where the elements are sorted
			% in the way expected by MPCController/construct()

			format = obj.xinitFormat;
			xinit = x0(:);

			% TODO: automatically assign default obj.xinitFormat when
			% importing from solvemp/polyunion
			if nargin==2 && (~isstruct(format) || length(format.names)==1)
				% only x0 required
				xinit = x0;
				return
			elseif mod(length(varargin), 2)~=0
				error('Arguments must come in key/value pairs.');
			end
			
			% parse varargin
			options = [];
			for i = 1:2:length(varargin)
				% validate each ['option', value] pair
				key = varargin{i};
				value = varargin{i+1};
				if ~ischar(key)
					error('Each option must be a string.');
				end
				% replace dots by underscores, e.g. "x.ref" -> "x_ref"
				key = strrep(key, '.', '_');
				% covert the key/value pair into a field of "options"
				options.(key) = value;
			end
			
			% now check that we have all required options
			%
			% obj.xinitFormat.names{1} always corresponds to 'x.init',
			% which we already have in xinit
			for i = 2:length(format.names)
				sane_name = strrep(format.names{i}, '.', '_');
				% validate presence
				if ~isfield(options, sane_name)
					% option is required
					error('Please provide initial value of "%s".', format.names{i});
				end
				value = options.(sane_name);
				% validate dimensions
				if ~isnumeric(value) || ~isequal(size(value), format.dims{i})
					error('"%s" must be a %dx%d vector.', format.names{i}, ...
						format.dims{i}(1), format.dims{i}(2));
				end
				% validation passed, include the initial value into xinit
				xinit = [xinit; value];
			end
		end
		
	end
	
end
