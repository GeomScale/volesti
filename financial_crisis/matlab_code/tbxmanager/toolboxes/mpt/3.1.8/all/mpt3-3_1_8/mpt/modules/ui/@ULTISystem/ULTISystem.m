classdef ULTISystem < LTISystem
    % Object representation of an uncertain LTI system
    %
    % To create an uncertain LTI system
    %   x^+ = A(L)x + B(L)u + E*d
    %     y = C(L)x + D(L)u
    % call:
    %   sys = ULTISystem('A', A, 'B', B, 'E', E, 'C', C, 'D', D)
    %
    % where "A", "B", "C", "D" can be cell arrays of vertex
    % representations of the uncertain dynamics.
    %
    % If "B", "C", "D" are not provided, they are
    % set to empty matrices. If "E" is not provided, E=eye(nx) is assumed.
    %
    % Note that the additive disturbance is only considered if the
    % "E" matrix is set not a non-zero and non-empty matrix.
    %
    % By default, the bounds on the disturbance are zero. These can be
    % changed via sys.d.min and sys.d.max.

    properties(SetAccess=protected)
        nd % number of additive disturbances
        E % Matrix multiplying disturbance in x(t+Ts) = A*x(t)+B*u(t)+E*d(t)
    end
	
    methods(Access = protected)

        function out = display_internal(obj)
            % Returns a string information about the system
            
            plural = @(s, n) [num2str(n) ' ' s repmat('s', 1, double(n~=1))];
            out = sprintf('%s with %s, %s, %s, %s', class(obj), ...
                plural('state', obj.nx), ...
                plural('input', obj.nu), ...
                plural('output', obj.ny), ...
                plural('disturbance', obj.nd));
        end

        function d = randomDisturbance(obj)
            % generates a random additive disturbance
            
            d = obj.d.boundsToPolyhedron().randomPoint();
        end
        
        function [rnd, lambdas] = randomConvexCombination(obj, middle)
            % generates a random convex combination of parametric
            % uncertainties
            
            if nargin==1
                % used to get a nominal model in the middle
                middle = false;
            end
            
            rnd = [];
            lambdas = [];
            
            function randomInstance(name)
                % random convex combination
                
                M = obj.(name);
                if ~iscell(M)
                    M = { M };
                end
                % random vector with 0<=L(:)<=1 and sum(L)=1
                if middle
                    L = ones(numel(M), 1)/numel(M);
                else
                    L = rand(numel(M));
                    L = L/sum(L);
                end
                Mr = zeros(size(M{1}));
                for i = 1:numel(M)
                    Mr = Mr + L(i)*M{i};
                end
                rnd.(name) = Mr;
                lambdas.(name) = L;
            end

            % these calls will update "rnd" and "lambdas"
            randomInstance('A');
            randomInstance('B');
            randomInstance('C');
            randomInstance('D');
        end
    end
    
	methods(Hidden)
		% implementation of abstract methods
		
		% no validation in these functions! it was already performed in
		% AbstractSystem/update() and output()

		function [xn, y, lambdas, d] = update_equation(obj, x, u)
			% returns the state update and the output

            % generate a random convex combination of parametric
            % uncertainties
            [rnd, lambdas] = obj.randomConvexCombination();
            d = [];
            
			xn = rnd.A*x;
            if obj.nd > 0
                d = obj.randomDisturbance();
                xn = xn + obj.E*d;
            end
            if ~isempty(u)
                xn = xn + rnd.B*u;
            end
			y = rnd.C*x;
            if ~isempty(u)
                y = y + rnd.D*u;
            end
		end

		function [y, lambdas] = output_equation(obj, x, u)
			% output equation
			
            % generate a random convex combination of parametric
            % uncertainties
            [rnd, lambdas] = obj.randomConvexCombination();
			y = rnd.C*x;
            if ~isempty(u)
                y = y + rnd.D*u;
            end
		end

		function out = has_feedthrough(obj)
			% feedthrough indication. must return true if the system has
			% direct feedthrough, false otherwise
			
            if iscell(obj.D)
                D = [obj.D{:}];
            else
                D = obj.D;
            end
			out = (nnz(D)~=0);
        end

        function sys = cellifyMatrices(obj)
            % Converts matrices to cells
            
            if ~iscell(obj.A)
                sys.A = { obj.A };
            else
                sys.A = obj.A;
            end
            if ~iscell(obj.B)
                sys.B = { obj.B };
            else
                sys.B = obj.B;
            end
            if ~iscell(obj.C)
                sys.C = { obj.C };
            else
                sys.C = obj.C;
            end
            if ~iscell(obj.D)
                sys.D = { obj.D };
            else
                sys.D = obj.D;
            end
        end

	end
	
    methods
        
        function obj = ULTISystem(varargin)
            % Constructor for uncertain LTI systems
            %
            % To create an uncertain LTI system
            %   x^+ = A(L)x + B(L)u + f + E*d
            %     y = C(L)x + D(L)u + g
            % call:
            %   s = ULTISystem('A', A, 'B', B, 'E', E, 'C', C, 'D', D, 'f', f, 'g', g)
            %
            % where "A", "B", "C", "D" can be cell arrays of vertex
            % representations of the uncertain dynamics.
            %
            % If "B", "C", "D", "f", "g" are not provided, they are
            % set to empty matrices.
            %
            % If "E" is not provided, E=eye(nx) is assumed.
            %
            % By default the disturbace "d" has zero lower/upper bounds.
            % These can be modified via sys.d.min and sys.d.max.
            
            if nargin == 0
                return
            end
            
            ip = inputParser;
            ip.KeepUnmatched = false;
            ip.addParamValue('A', []);
            ip.addParamValue('B', []);
            ip.addParamValue('f', [], @validate_realmatrix);
            ip.addParamValue('E', [], @validate_realmatrix);
            ip.addParamValue('C', []);
            ip.addParamValue('D', []);
            ip.addParamValue('g', [], @validate_realmatrix);
            ip.addParamValue('domain', [], @validate_polyhedron);
            ip.addParamValue('Ts', 1, @(x) isa(x, 'double') && numel(x)==1 && x>=0);
            ip.parse(varargin{:});
            S = ip.Results;
                
            obj.A = S.A;
            obj.B = S.B;
            obj.E = S.E;
            obj.C = S.C;
            obj.D = S.D;
            try
                obj.Ts = S.Ts;
            end
            
            if iscell(obj.A)
                obj.nx = size(obj.A{1}, 2);
            else
                obj.nx = size(obj.A, 2);
            end
            if iscell(obj.B)
                obj.nu = size(obj.B{1}, 2);
            else
                obj.nu = size(obj.B, 2);
            end
            if iscell(obj.C)
                obj.ny = size(obj.C{1}, 1);
            else
                obj.ny = size(obj.C, 1);
            end
            if isempty(obj.E)
                obj.E = eye(obj.nx);
            end
            obj.nd = size(obj.E, 2);

			if isempty(obj.B)
				obj.B = zeros(obj.nx, obj.nu);
            end
            if isempty(obj.E)
                obj.E = zeros(obj.nx, obj.nd);
            end
			if isempty(obj.C)
				obj.C = zeros(obj.ny, obj.nx);
			end
			if isempty(obj.D)
				obj.D = zeros(obj.ny, obj.nu);
            end
            
            % TODO: validate dimensions of state-update/output matrices
			if size(obj.E, 1)~=obj.nx
                error('"E" must have %d row(s).', obj.nx);
            end

            % TODO: unify adding of signals/domains with LTISystem
            x = SystemSignal(obj.nx);
            x.name = 'x';
            x.setKind('x');
            obj.addComponent('x', x);
            
            u = SystemSignal(obj.nu);
            u.name = 'u';
            u.setKind('u');
            obj.addComponent('u', u);
            
            y = SystemSignal(obj.ny);
            y.name = 'y';
            y.setKind('y');
            obj.addComponent('y', y);
            
            if obj.nd > 0
                d = SystemSignal(obj.nd);
                d.without('penalty');
                d.name = 'd';
                d.setKind('d');
                obj.addComponent('d', d);
                obj.d.min = zeros(obj.d.n, 1);
                obj.d.max = zeros(obj.d.n, 1);
            end
               
			% domain is a polyhedron in the x-u space
			if isobject(S.domain)
				obj.domain = S.domain;
			else
				obj.domain = obj.x.boundsToPolyhedron()*obj.u.boundsToPolyhedron();
			end
			
			% domain in the x-space
			
			% TODO: remove the hard-coded 'fourier' method once
			% Polyhedron/projection is reliable
			obj.domainx = obj.domain.projection(1:obj.nx, 'fourier');

			if nargin==1 && isstruct(varargin{1})
                % import constraints from sysStruct
                obj.importSysStructConstraints(S);
			end

        end
        
        function out = toSS(obj)
            % Converts the ULTI system into a state-space object
            error('Conversion to a state-space object not supported.');
        end
        
		function [K, P, dyn] = LQRGain(obj)
            % Computes the LQR gain
            %
            % For an uncertain system x^+ = A(L)*x + B(L)*u this method
            % computes the gain K such that the closed-loop system
            % x^+ = (A(L)+B(L)*K)*x is stable for all realizations of the
            % parametric uncertainty.
            %
            % Additive disturbances are not supported.
            
            % system must have inputs
            assert(obj.nu>0, 'The system must have control inputs.');
            % penalties must be defined
            assert(obj.x.hasFilter('penalty') && isa(obj.x.penalty, 'QuadFunction'), ...
                'The state penalty must be a quadratic function.');
            assert(obj.u.hasFilter('penalty') && isa(obj.u.penalty, 'QuadFunction'), ...
                'The input penalty must be a quadratic function.');
            
            [K, P] = obj.stabilizingGain(obj.x.penalty.H, obj.u.penalty.H);
        end
        
        function [K, P] = stabilizingGain(obj, Qx, Qu)
            % Computes the stabilizing state-feedback gain
            %
            % For an uncertain system x^+ = A(L)*x + B(L)*u this method
            % computes the gain K such that the closed-loop system
            % x^+ = (A(L)+B(L)*K)*x is stable for all realizations of the
            % parametric uncertainty.
            %
            % Additive disturbances are not supported.
            
            global MPTOPTIONS
            if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
            
            % system must have inputs
            assert(obj.nu>0, 'The system must have control inputs.');

            % this method also computes the LQR gain if input/state
            % penalties are defined
            if nargin<=1
                Qx = [];
                Qu = [];
            end
            
            % the disturbance set must be either empty or just zero
%             D = obj.d.boundsToPolyhedron();
%             assert(D.isEmptySet() || (~D.isFullDim() && D.contains(zeros(obj.nd, 1))), ...
%                 'Additive disturbances are not supported.');
            
            % convert all matrices to cells
            unc = obj.cellifyMatrices();
            
            % set up the LMI
            %
            % also see: http://robot2.disp.uniroma2.it/~zack/ftp/users/DeSimone/dilated1.pdf
            Q = sdpvar(obj.nx, obj.nx, 'symmetric');
            Y = sdpvar(obj.nu, obj.nx, 'full');
            prob = [0 <= Q <= MPTOPTIONS.infbound];
            for ia = 1:numel(unc.A)
                for ib = 1:numel(unc.B)
                    if isempty(Qx)
                        prob = prob + [ [Q, (unc.A{ia}*Q+unc.B{ib}*Y)'; ...
                            (unc.A{ia}*Q+unc.B{ib}*Y), Q ] >= 0 ];
                    else
                        prob = prob + [ [Q, (unc.A{ia}*Q+unc.B{ib}*Y)', (Qx^0.5*Q)', (Qu^0.5*Y)'; ...
                            (unc.A{ia}*Q+unc.B{ib}*Y), Q, zeros(obj.nx), zeros(obj.nx, obj.nu); ...
                            (Qx^0.5*Q), zeros(obj.nx), eye(obj.nx), zeros(obj.nx, obj.nu); ...
                            (Qu^0.5*Y), zeros(obj.nu, obj.nx), zeros(obj.nu, obj.nx), eye(obj.nu)] >= 0 ];
                    end
                end
            end
            options = sdpsettings('verbose', MPTOPTIONS.verbose);
            sol = solvesdp(prob, -trace(Q), options);
            % these statuses indicate a feasible solution
            % (from "help yalmiperror"):
            %  0: Successfully solved
            %  3: Maximum #iterations or time-limit exceeded
            %  4: Numerical problems
            %  5: Lack of progress
            if ismember(sol.problem, [3, 4, 5])
                warning(sol.info);
            elseif sol.problem~=0
                error(sol.info);
            end
            
            P = inv(double(Q));
            K = double(Y)*P;
        end
		
		function P = LQRPenalty(obj)
			% Returns the LQR penalty
            
            [~, P] = obj.LQRGain();
		end
		
		function S = LQRSet(obj)
			% Returns the LQR invariant set
            
            K = obj.LQRGain();
            loop = K*obj;
            S = loop.toSystem().invariantSet();
        end
		
        function C = LQRController(obj)
			% Computes a stabilizing LQR controller

            [K, P] = obj.LQRGain();
            L = LQRController(obj, K, P).toInvariant();
            C = EMPCController(L);
        end
        
        function C = stabilizingController(obj)
			% Computes a stabilizing controller

            [K, P] = obj.stabilizingGain();
            L = SFController(obj, K, P).toInvariant();
            C = EMPCController(L);
        end
		
		function Z = reachableSet(obj, varargin)
			% Computes the forward/backwards reachable sets
            %
            % Consider a system x^+ = g(x, u, d) = A(L)*x + B(L)*u + E*d 
            % where A(L) = \sum L_i*A_i, 0<=L_i<=1, \sum L_i = 1 and B(L)
            % defined similarly. Then
            %
            % sys.reachableSet('direction', 'backward', 'X', X, 'U', U, 'D', D)
            % computes the one-step pre-set of a X, i.e.
            %
            % Pre(X) = { x | \exists u \in U s.t. g(x, u, d) \in X, \forall d \in D }
            %
            % If X, U, and/or D are omitted, the sets are created based on
            % lower/upper bounds in sys.x.min/max, sys.u.min/max and
            % sys.d.min/max, respectively.
            %
            % The method also supports autonomous systems with no control
            % inputs.
            %
            % sys.reachableSet('direction', 'forward', 'X', X, 'U', U, 'D', D)
            % computes the one-step reach set of the set X, i.e.,
            %
            % Reach(X) = { g(x, u, d) | \exists x \in X, u \in U, d \in D }
            %
            % Again, the X, U, and/or D inputs can be omitted.
            %
            % If the "direction" argument is omitted, "forward" is the
            % default. Forward reach sets can also be computed for
            % autonomous systems.
            %
            % WARNING: If the dynamics is subject to parametric
            % uncertainties, then the forward reachable set computed by
            % this method is just a coarse approximation of the true set
            % computed for vertices of the parametric uncertainty. A more
            % reasonable inner approximation can be obtained by manually
            % specifying a fine grid of the parametric uncertainty using
            % the "A" and "B" optional inputs.

            ip = inputParser;
            ip.KeepUnmatched = false;
            ip.addParamValue('X', [], @validate_polyhedron);
            ip.addParamValue('U', [], @validate_polyhedron);
            ip.addParamValue('D', [], @validate_polyhedron);
            ip.addParamValue('direction', 'forward', @ischar);
            ip.addParamValue('A', {}, @iscell);
            ip.addParamValue('B', {}, @iscell);
            ip.parse(varargin{:});
            options = ip.Results;
            
            if isempty(options.X)
                options.X = obj.x.boundsToPolyhedron() & obj.domainx;
            end
            if isempty(options.U)
                options.U = obj.u.boundsToPolyhedron();
            end
            if isempty(options.D)
                options.D = obj.d.boundsToPolyhedron();
            end
            
            assert(numel(options.X)==1, '"X" must be a single polyhedron.');
            assert(options.X.Dim==obj.nx, sprintf('"X" must be in %dD.', obj.nx));
            assert(numel(options.D)==1, '"D" must be a single polyhedron.');
            assert(options.D.Dim==obj.nd, sprintf('"D" must be in %dD.', obj.nd));
            assert(numel(options.U)==1, '"U" must be a single polyhedron.');
            assert(options.U.Dim==obj.nu, sprintf('"U" must be in %dD.', obj.nu));

            if options.X.isEmptySet()
                % exit quickly if the target is empty
                Z = options.X.copy();
                return
            end

            backwards = lower(options.direction(1))=='b';
            param_unc = iscell(obj.A) || iscell(obj.B);

            function Z = preset(A, B)
                % pre set (backward)
                
                if obj.nu==0
                    % autonomous systems
                    Z = invAffineMap(options.X-(obj.E*options.D), A);
                else
                    % non-autonomous systems
                    Z = invAffineMap((options.X-(obj.E*options.D))+((-B)*options.U), A);
                end
            end
            
            function Z = reachset(A, B)
                % reach set (forward)
                
                if obj.nu==0
                    % autonomous systems
                    Z = (A*options.X)+(obj.E*options.D);
                else
                    Z = (A*options.X)+(B*options.U)+(obj.E*options.D);
                end
            end

            % convert all to cells to simply the code
            if isempty(options.A)
                A = obj.A;
            else
                % use user-specified grid of the uncertainty polytope
                A = options.A;
            end
            if ~iscell(A)
                A = { A };
            end
            if isempty(options.B)
                B = obj.B;
            else
                % use user-specified grid of the uncertainty polytope
                B = options.B;
            end
            if ~iscell(B)
                B = { B };
            end
            
            if backwards
                % preset as the intersection of all presets
                % important: must project the final intersection, not
                % intermediate results
                H = []; K = [];
                D = (options.X-(obj.E*options.D));
                if obj.nu>0
                    D = D*options.U;
                end
                for i = 1:numel(A)
                    for j = 1:numel(B)
                        if obj.nu>0
                            Acl = [A{i}, B{j}; zeros(obj.nu, obj.nx), eye(obj.nu)];
                        else
                            Acl = A{i};
                        end
                        H = [H; D.A*Acl];
                        K = [K; D.b];
                    end
                end
                Z = Polyhedron(H, K);
                if obj.nu>0
                    Z = Z.projection(1:obj.nx);
                end
            else
                Z = reachset(A{1}, B{1});
                % union of all combinations
                for i = 1:numel(A)
                    for j = 1:numel(B)
                        if i==1 && j==1
                            % this was already covered above
                            continue
                        end
                        Q = reachset(A{i}, B{j});
                        % reach set is the union
                        if ~Q.isEmptySet()
                            Z = [Z, Q];
                        end
                    end
                end
            end
		end

		function [Z, converged, iters] = invariantSet(obj, varargin)
			% Computes invariant set of the system
            %
            % Consider an autonomous system x(k+1) = A(L)*x(k) + E*d(k) with
            % x \in X and d \in D.
            %
            % Then sys.invariantSet() computes the maximal robust positive
            % invariant set Oinf defined by
            %   Oinf = { x(0) \in X | x(k) \in Oinf, \forall d \in D, k>0 }
            %
            % For the system x(k+1) = A(L)*x(k) + B(L)*u(k) + E*d(k) with
            % x \in X, u \in U and d \in D, sys.invariantSet() computes the
            % maximal robust control invariant set Cinf
            %   Cinf = { x(0) \in X | \exists u(k) \in U s.t. x(k+1) \in
            %                         Cinf, \forall d \in D, \forall k>0 }
            %
            % The long syntax is
            %  sys.invariantSet('X', X, 'U', U, 'D', D, 'maxIterations', m)
            % where "X", "U" and "D" are polytopes bounding the respective
            % quantities. If omitted, these bounds are extracted from
            % min/max bounds of respective signals (sys.x, sys.u and sys.d,
            % respectively).
            %
            % The "maxIterations" input allows to prematurely abort the
            % recursive computation. Note that in such a case the returned
            % set is not necessarily maximal.

			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
           
			ip = inputParser;
			ip.KeepUnmatched = false;
            ip.addParamValue('maxIterations', ...
				MPTOPTIONS.modules.ui.invariantSet.maxIterations, ...
                @isscalar);
            ip.addParamValue('U', [], @validate_polyhedron);
            ip.addParamValue('X', [], @validate_polyhedron);
            ip.addParamValue('D', [], @validate_polyhedron);
			ip.parse(varargin{:});
			options = ip.Results;

            if isempty(options.X)
                options.X = obj.x.boundsToPolyhedron() & obj.domainx;
            end
            if isempty(options.U)
                options.U = obj.u.boundsToPolyhedron();
            end
            if isempty(options.D)
                options.D = obj.d.boundsToPolyhedron();
            end

            assert(numel(options.D)==1, '"D" must be a single polyhedron.');
            assert(options.D.Dim==obj.nd, sprintf('"D" must be in %dD.', obj.nd));
            assert(numel(options.U)==1, '"U" must be a single polyhedron.');
            assert(options.U.Dim==obj.nu, sprintf('"U" must be in %dD.', obj.nu));
            assert(numel(options.X)==1, '"X" must be a single polyhedron.');
            assert(options.X.Dim==obj.nx, sprintf('"X" must be in %dD', obj.nx));
            
            Zo = options.X;
			converged = false;
            for i = 1:options.maxIterations
                fprintf('Iteration %d...\n', i);
                Zn = obj.reachableSet('X', Zo, 'direction', 'backward', ...
                    'U', options.U, 'D', options.D) & Zo;
                if Zn==Zo
                    converged = true;
                    break
                else
                    Zo = Zn;
                end
            end
            iters = i;
            if ~converged
                warning('Computation finished without convergence.');
            end
            Z = Zn;
		end
                
		function con = constraints(obj)
			% Convert LTI model into YALMIP constraints

			% constraints on variables
			con = constraints@AbstractSystem(obj);

            function add_x_constraint(a_xn, a_xp, a_u, a_A, a_B)
                con = con + [ a_xn == a_A*a_xp + a_B*a_u ];
            end
            function add_y_constraint(a_y, a_xp, a_u, a_C, a_D)
                con = con + [ a_y == a_C*a_xp + a_D*a_u ];
            end

			% add the LTI dynamics constraints
			x = obj.x.var;
			u = obj.u.var;
			y = obj.y.var;
            
            % robustify state constraints
            Xrob = obj.x.boundsToPolyhedron() - obj.E*obj.d.boundsToPolyhedron();

            % do we have parametric uncertainty?
            param_unc = iscell(obj.A) || iscell(obj.B) || ...
                iscell(obj.C) || iscell(obj.D);
                
            if param_unc
                % get the nominal model
                nominal = obj.randomConvexCombination(true);
                % convert all system matrices to cells
                unc = obj.cellifyMatrices();
                
                % first set up the problem for nominal dynamics (this will
                % be used to set up the cost function as well)
                for k = 1:obj.Internal.system.N
                    if obj.nx > 0
                        add_x_constraint(x(:, k+1), x(:, k), u(:, k), ...
                            nominal.A, nominal.B);
                    end
                    if obj.ny > 0
                        add_y_constraint(y(:, k), x(:, k), u(:, k), ...
                            nominal.C, nominal.D);
                    end
                    % robust state constraints
                    con = con + [ ismember(x(:, k), Xrob) ];
                end
                
                % now introduce additional variables for each vertex of the
                % uncertainty polytope
                xprevious = { x(:, 1) };
                if obj.ny > 0
                    Yb = obj.y.boundsToPolyhedron();
                end
                
                for k = 1:obj.Internal.system.N
                    xnew = {};
                    for j = 1:numel(xprevious)
                        xp = xprevious{j};
                        con = con + [ ismember(xp, Xrob) ];
                        if obj.nx > 0
                            for ia = 1:numel(unc.A)
                                for ib = 1:numel(unc.B)
                                    xn = sdpvar(obj.nx, 1);
                                    add_x_constraint(xn, xp, u(:, k), ...
                                        unc.A{ia}, unc.B{ib});
                                    xnew{end+1} = xn;
                                end
                            end
                        end
                        if obj.ny > 0
                            for ic = 1:numel(unc.C)
                                for id = 1:numel(unc.D)
                                    yn = sdpvar(obj.ny, 1);
                                    add_y_constraint(yn, xp, u(:, k), ...
                                        unc.C{ic}, unc.D{id});
                                    con = con + [ ismember(yn, Yb) ];
                                end
                            end
                        end
                    end
                    xprevious = xnew;
                end

                % add terminal set constraints
                if obj.x.hasFilter('terminalSet')
                    for i = 1:numel(xprevious)
                        % at this point xprevious contains all possible
                        % realizations of x(:, N+1)
                        con = con + [ ismember(xprevious{i}, obj.x.terminalSet) ];
                    end
                end
                
            else
                % no parametric uncertainties, just use tightened state
                % constraints
                for k = 1:obj.Internal.system.N
                    if obj.nx > 0
                        add_x_constraint(x(:, k+1), x(:, k), u(:, k), ...
                            obj.A, obj.B);
                    end
                    if obj.ny > 0
                        add_y_constraint(y(:, k), x(:, k), u(:, k), ...
                            obj.C, obj.D);
                    end
                    % robust state constraints
                    con = con + [ ismember(x(:, k), Xrob) ];
                end
            end
		end
	end

end
