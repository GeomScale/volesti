classdef PWASystem < AbstractSystem
%
%  PWASYSTEM: Represents discrete-time piecewise affine systems 
%  =============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      sys = PWASystem([lti1, lti2, ..., ltiM])
%    
%  
%  DESCRIPTION
%  -----------
%     This class represents PWA systems, which are composed of a finite number of
%  local affine dynamics, each valid in a corresponding polyhedral region of the
%  state-input space: 
%                                                                            
%                x(t+Ts)  = A  x(t) + B  u(t) + f   if  (x,u) in R           
%                            i         i         i                i          
%                                                                            
%                   y(t)  = C  x(t) + D  u(t) + g                            
%                            i         i         i                           
%     where x in R^n_x  is the state vector, u in R^n_u  is the vector of inputs, y
%  in R^n_y  is the vector of outputs, T_srepresents the sampling time, and R_i 
%  are the polyhedral regions of validity of the i-th local dynamics.
%    Each PWA system defines following properties: 
%    
%     - A, B, f: matrices of the state-update equation, stored as cell arrays
%     (read-only) 
%     - C, D, g: matrices of the output equation, stored as cell arrays (read-only)
%     
%     - Ts: sampling time (read-only) 
%     - domain: array of polyhedra denoting domain of the i-th local model
%     (read-only) 
%     - nx, nu, ny: number of states, inputs and outputs (automatically determined,
%     read-only) 
%     - ndyn: number of local models (read-only) 
%     - x: specifications of system's states (see help SystemSignal) 
%     - u: specifications of system's inputs (see help SystemSignal) 
%     - y: specifications of system's outputs (see help SystemSignal) 
%     - d: specifications of the binary dynamics selector signal (see help
%     SystemSignal) 
%    The preferred way to define a PWA system consisting of a finite number of
%  local affine models is to provide the list of LTI models to the PWASystem
%  constructor:
%     pwasys = PWASystem([ltisys1, ltisys2, ..., ltisysM])
%    Here, each LTI model must have its domain defined by the ltisys.setDomain()
%  method (see " help LTISystem/setDomain").
%  
%  SEE ALSO
%  --------
%     LTISystem
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
 
 
    properties(SetAccess=protected)
        A % Matrices of the state-update equation x(t+Ts)=A_i*x(t)+B_i*u(t)+f_i
        B % Matrices of the state-update equation x(t+Ts)=A_i*x(t)+B_i*u(t)+f_i
        C % Matrices of the output equation y(t)=C_i*x(t)+D_i*u(t)+g_i
        D % Matrices of the output equation y(t)=C_i*x(t)+D_i*u(t)+g_i
        f % Matrices of the state-update equation x(t+Ts)=A_i*x(t)+B_i*u(t)+f_i
        g % Matrices of the output equation y(t)=C_i*x(t)+D_i*u(t)+g_i
        ndyn % Number of local dynamics
	end

	properties(SetAccess=private, GetAccess=private, Hidden)
		modes % local affine modes as an array of LTISystem
		polyunion % state-update and output equations as a PolyUnion object
	end
            
	methods(Hidden)
		% implementation of abstract methods
		
		% no validation in these functions! it was already performed in
		% AbstractSystem/update() and output()
		
		function [xn, y] = update_equation(obj, x, u)
			% returns the state update and the output
			
			% use the first-region tiebreak
			xn = obj.polyunion.feval([x; u], 'update', 'tiebreak', @(q) 0);
			y = obj.polyunion.feval([x; u], 'output', 'tiebreak', @(q) 0);
		end
		
		function y = output_equation(obj, x, u)
			% output equation

			% use the first-region tiebreak
			y = obj.polyunion.feval([x; u], 'output', 'tiebreak', @(q) 0);
		end
		
		function out = has_feedthrough(obj)
			% feedthrough indication. must return true if the system has
			% direct feedthrough, false otherwise

			out = (nnz(cat(2, obj.D{:}))~=0);
		end
    end
    
    methods(Access=protected)
        
        function out = display_internal(obj)
            % Returns a string information about the system
            
            plural = @(s, n) [num2str(n) ' ' s repmat('s', 1, double(n~=1))];
            out = sprintf('%s with %s, %s, %s, %s', class(obj), ...
                plural('state', obj.nx), ...
                plural('input', obj.nu), ...
                plural('output', obj.ny), ...
                plural('mode', obj.ndyn));
        end
    end


    methods
        
        function obj = PWASystem(varargin)
            % Constructor for PWA systems
            %
            %   x^+ = A{i}*x + B{i}*u + f{i}   if   [x; u] \in P(i)
            %     y = C{i}*x + D{i}*u + g{i}
            % 
            % s = PWASystem('A', A, 'B', B, 'C', C, 'D', D, 'domain', P)
			%
			% To import from a sysStruct structure:
            %   s = PWASystem(sysStruct)

            
            if nargin == 0
                return
			end
            
			if nargin==1 && isstruct(varargin{1})
                % import from sysStruct
                S = mpt_verifySysStruct(varargin{1});
                P = [];
                for i = 1:length(S.guardX)
                    P = [P, Polyhedron('A', [S.guardX{i}, S.guardU{i}], ...
                        'b', S.guardC{i}).minHRep];
                end
                S.domain = P;
				
			elseif nargin==1 && isa(varargin{1}, 'MLDSystem')
				% convert MLD to PWA
				M = varargin{1};
				obj = M.toPWA();
				return
				
			elseif nargin==1 && isa(varargin{1}, 'LTISystem')
				% convert array of LTI systems into PWA
				
				lti = varargin{1};
				S = [];
				S.A = cell(1, numel(lti));
				S.B = cell(1, numel(lti));
				S.f = cell(1, numel(lti));
				S.C = cell(1, numel(lti));
				S.D = cell(1, numel(lti));
				S.g = cell(1, numel(lti));
				
				% check that all systems have identical dimensions
				nx = cell(1, numel(lti));
				nu = cell(1, numel(lti));
				ny = cell(1, numel(lti));
				[nx{:}] = lti.nx; nx = cat(2, nx{:});
				[nu{:}] = lti.nu; nu = cat(2, nu{:});
				[ny{:}] = lti.ny; ny = cat(2, ny{:});
				if any(diff(nx))
					error('All systems must have identical state dimensions.');
				end
				if any(diff(nu))
					error('All systems must have identical input dimensions.');
				end
				if any(diff(ny))
					error('All systems must have identical output dimensions.');
				end
				
				% check that all systems have identical sampling time
				S.Ts = lti(1).Ts;
				for i = 2:numel(lti)
					if ~isequal(lti(i).Ts, S.Ts)
						error('All systems must have identical sampling time.');
					end
				end

				% put matrices into a cell array
				for i = 1:numel(lti)
					S.A{i} = lti(i).A;
					S.B{i} = lti(i).B;
					S.f{i} = lti(i).f;
					S.C{i} = lti(i).C;
					S.D{i} = lti(i).D;
					S.g{i} = lti(i).g;
				end
				
				% warn user that constraints will be ignored
				for i = 1:numel(lti)
					if hasFilter(lti(i).x, 'min') || ...
							hasFilter(lti(i).x, 'max') || ...
							hasFilter(lti(i).u, 'min') || ...
							hasFilter(lti(i).u, 'max') || ...
							hasFilter(lti(i).y, 'min') || ...
							hasFilter(lti(i).y, 'max')
						% TODO: set constraints automatically
						fprintf('State/input/output constraints not imported, set them manually afterwards.\n');
						break;
					end
				end
				
				% import domain
				S.domain = [];
				for i = 1:numel(lti)
					S.domain = [S.domain; lti(i).domain];
				end
				
			elseif nargin==1
				error('Unrecognized import source.');
			
			else
				% import from option/value pairs
				ip = inputParser;
				ip.KeepUnmatched = false;
				ip.addParamValue('A', [], @iscell);
				ip.addParamValue('B', [], @iscell);
				ip.addParamValue('f', [], @iscell);
				ip.addParamValue('C', [], @iscell);
				ip.addParamValue('D', [], @iscell);
				ip.addParamValue('g', [], @iscell);
				ip.addParamValue('domain', [], @validate_polyhedron);
				ip.addParamValue('Ts', 1, @(x) isa(x, 'double') && numel(x)==1 && x>=0);
				ip.addParamValue('sysStruct', [], @isstruct);
				ip.parse(varargin{:});
				S = ip.Results;
			end

            obj.A = S.A;
            obj.B = S.B;
            obj.f = S.f;
            obj.C = S.C;
            obj.D = S.D;
            obj.g = S.g;
            obj.Ts = S.Ts;
			
			% TODO: check that domains in the x-u space do not overlap
            obj.domain = S.domain;

            obj.nx = size(obj.A{1}, 2);
            obj.nu = size(obj.B{1}, 2);
            obj.ny = size(obj.C{1}, 1);
            obj.ndyn = length(obj.A);

			obj.domainx = [];
			for i = 1:length(obj.domain)
				% TODO: remove the hard-coded 'fourier' method once
				% Polyhedron/projection is reliable
				dx = obj.domain(i).projection(1:obj.nx, 'fourier').minHRep();
				obj.domainx = [obj.domainx dx];
			end
            
            if isempty(obj.f)
                obj.f = zeros(obj.nx, 0);
            end
            if isempty(obj.g)
                obj.g = zeros(obj.nx, 0);
            end

            x = SystemSignal(obj.nx);
            x.name = 'x';
            x.setKind('x');
            obj.addComponent('x', x);
            
            u = SystemSignal(obj.nu);
            u.name = 'u';
            u.setKind('u');
            obj.addComponent('u', u);
			if isfield(S, 'Uset')
				% mark inputs as binary
				obj.u.with('binary');
				obj.u.binary = find(cellfun(@(x) isequal(x, [0 1]), S.Uset));
			end

			y = SystemSignal(obj.ny);
			y.name = 'y';
			y.setKind('y')
			obj.addComponent('y', y);
            
            % create additional binary selectors for the PWA dynamics
            d = SystemSignal(obj.ndyn);
			d.addFilter('binary');
			d.binary = true;
            d.name = 'd';
            d.setKind('d');
            obj.addComponent('d', d);
            
            if nargin==1
                % import constraints from sysStruct
                obj.importSysStructConstraints(S);
			end
			
			% store the state-update and output functions as a PolyUnion
			P = Polyhedron(obj.domain); % create copy
			for i = 1:obj.ndyn
				update = AffFunction([obj.A{i}, obj.B{i}], obj.f{i});
				output = AffFunction([obj.C{i}, obj.D{i}], obj.g{i});
				P(i).addFunction(update, 'update');
				P(i).addFunction(output, 'output');
			end
			obj.polyunion = PolyUnion(P);
			
			% store the individual modes as an array of LTISystems
			modes = [];
			for pwa_index = 1:obj.ndyn
				lti = LTISystem('A', obj.A{pwa_index}, 'B', obj.B{pwa_index}, ...
					'C', obj.C{pwa_index}, 'D', obj.D{pwa_index}, ...
					'f', obj.f{pwa_index}, 'g', obj.g{pwa_index}, ...
					'domain', obj.domain(pwa_index), 'Ts', obj.Ts);
				lti.x = obj.x;
				lti.u = obj.u;
				lti.y = obj.y;
				modes = [modes, lti];
			end
			obj.modes = modes;
		end
        
		function out = toLTI(obj, pwa_index)
			% Converts dynamics indexed by "pwa_index" to an LTI system
			
			if pwa_index<1 || pwa_index>obj.ndyn
				error('Index out of range.');
			end
			out = obj.modes(pwa_index);
		end
		
		function [S, SN, dyn, dynN] = reachableSet(obj, varargin)
			% Computes the forward/backwards reachable N-step set

			ip = inputParser;
			ip.KeepUnmatched = false;
			ip.addParamValue('direction', 'forward', @ischar);
			ip.addParamValue('N', 1, @isnumeric);
            ip.addParamValue('X', ...
				[], ...
				@validate_polyhedron);
            ip.addParamValue('U', ...
				[], ...
				@validate_polyhedron);
			ip.addParamValue('merge', true, @islogical);
			ip.parse(varargin{:});
			options = ip.Results;

            if isempty(options.X)
                options.X = obj.x.boundsToPolyhedron();
            end
            if isempty(options.U)
                options.U = obj.u.boundsToPolyhedron();
            end
            
			if numel(options.U)~=1
				error('Input constraints must be a single polyhedron.');
			end
			if obj.nu>0 && options.U.isEmptySet()
				error('Input constraints must not be empty.');
			end
			if any(arrayfun(@(x) x.Dim~=obj.nx, options.X))
				error('State constraints must be a polyhedron in %dD.', obj.nx);
			end
			if obj.nu>0 && options.U.Dim~=obj.nu
				error('Input constraints must be a polyhedron in %dD.', obj.nu);
            end

            % only keep non-empty sets
            X = options.X(find(~[options.X.isEmptySet()]));
            if isempty(X)
                % all targets are empty -> reach set is empty and exit
                S = Polyhedron.emptySet(obj.nx);
                SN{1} = S;
                dyn = 0;
                return
            end
			SN = {}; dynN = {0};
			for n = 1:options.N
                Xp = []; dyn = [];
				for j = 1:obj.ndyn
					lti = obj.toLTI(j);
					R = lti.reachableSet('N', 1, ...
						'direction', options.direction, ...
						'X', X, 'U', options.U);
					% is the union of "R" convex?
					if numel(R)>1
						H = PolyUnion(R);
						if H.isConvex
							R = H.convexHull;
						elseif options.merge
							R = PolyUnion(R).merge().Set;
						end
					end
					Xp = [Xp, R];
					dyn = [dyn j*ones(1, numel(R))];
                end
                keep = find(~Xp.isEmptySet());
                SN{n} = Xp(keep);
                dynN{n} = dyn(keep);
                X = SN{end};
                if isempty(X)
                    break
                end
            end
            if isempty(SN)
                SN{1} = Polyhedron.emptySet(obj.nx);
            end
			S = SN{end};
            % only keep non-empty sets
            S = S(find(~[S.isEmptySet()]));
            if isempty(S)
                S = Polyhedron.emptySet(obj.nx);
            end
			dyn = dynN{end};
		end

		function [X, dyn] = invariantSet(obj, varargin)
			% Computes invariant set of the system
			

			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			options = MPTOPTIONS.modules.ui.invariantSet;

			ip = inputParser;
			ip.KeepUnmatched = false;
            ip.addParamValue('maxIterations', ...
				options.maxIterations, @isscalar);
            ip.addParamValue('X', ...
				obj.domainx, ...
				@validate_polyhedron);
            ip.addParamValue('U', ...
				obj.u.boundsToPolyhedron(), ...
				@validate_polyhedron);
			ip.addParamValue('merge', true, @islogical);
			ip.parse(varargin{:});
			options = ip.Results;

			if numel(options.U)~=1
				error('Input constraints must be a single polyhedron.');
			end
			if obj.nu>0 && options.U.isEmptySet()
				error('Input constraints must not be empty.');
			end
			if any(arrayfun(@(x) x.Dim~=obj.nx, options.X))
				error('State constraints must be a polyhedron in %dD.', obj.nx);
			end
			if obj.nu>0 && options.U.Dim~=obj.nu
				error('Input constraints must be a polyhedron in %dD.', obj.nu);
			end

			Xo = options.X;
			U = options.U;
			converged = false;
			for i = 1:options.maxIterations
				fprintf('Iteration %d...\n', i);
				[X, ~, dyn] = obj.reachableSet('X', Xo, 'U', U, ...
					'direction', 'backward', ...
					'merge', true);
				X = X.intersect(obj.x.boundsToPolyhedron);
				if all(X.isEmptySet()) || X==Xo
					converged = true;
					break
				else
					Xo = X;
				end
			end
			if ~converged
				warning('Computation finished without convergence.');
			end
		end

        function C = constraints(obj, k)
            % Convert LTI model into YALMIP constraints
            
            % constraints on variables
            C = constraints@AbstractSystem(obj);
            % add the PWA dynamics constraints
            x = obj.x.var;
            u = obj.u.var;
            y = obj.y.var;
            d = obj.d.var;
            for k = 1:obj.Internal.system.N
                if obj.nx > 0
                    for dyn = 1:obj.ndyn
                        C = C + [ implies(d(dyn, k), ...
                            x(:, k+1) == obj.A{dyn}*x(:, k) + obj.B{dyn}*u(:, k) + obj.f{dyn}) ];
                    end
                end
                if obj.ny > 0
                    for dyn = 1:obj.ndyn
                        C = C + [ implies(d(dyn, k), ...
                            y(:, k) == obj.C{dyn}*x(:, k) + obj.D{dyn}*u(:, k) + obj.g{dyn}) ];
                    end
                end
                for dyn = 1:obj.ndyn
					if isa(obj.domain(dyn), 'polytope')
						[H, K] = double(obj.domain(dyn));
					else
						H = obj.domain(dyn).A;
						K = obj.domain(dyn).b;
						H = [H; obj.domain(dyn).Ae; -obj.domain(dyn).Ae];
						K = [K; obj.domain(dyn).be; -obj.domain(dyn).be];
					end
                    C = C + [ implies(d(dyn, k), H*[x(:, k); u(:, k)] <= K) ];
                end
                C = C + [ sum(d(:, k)) == 1 ];
            end
        end

        function map = transitionMap(obj)
            % Transition map for an autonomous PWA system
            %
            % Given a PWA system x^+ = A_i*x+f_i IF x \in R_i, this method
            % computes the transition map M as an n-by-n logical matrix
            % with M(i, j)=true iff \exists x \in \R_i such that
            % x^+ \in \R_j. Also returns the subset of \R_i which enters
            % \R_j.
            %
            % Limitation: only full-dimensional transitions are considered
            %
            % Syntax:
            %
            %   map = pwa.transitionMap()
            %
            % Input:
            %   pwa: autonomous PWASystem object with
            %           x^+ = A_i*x+f_i IF x \in R_i
            %
            % Output:
            %   map.transition: n-by-n logical matrix with the i-j entry
            %                   set to true iff there exists a transition
            %                   from the i-th region to the j-th region
            %      map.regions: n-by-n cell array of polyhedra with the i-j
            %                   entry representing the subset of R_i which
            %                   enters R_j

            global MPTOPTIONS

            error(obj.rejectArray());
            if obj.nu>0
                error('This method only supports autonomous systems.');
            end

            map.transitions = false(obj.ndyn, obj.ndyn);
            map.regions = cell(obj.ndyn, obj.ndyn);
            tic;
            for i = 1:obj.ndyn
                if toc > MPTOPTIONS.report_period
                    fprintf('progress: %d/%d\n', i, obj.ndyn);
                    tic;
                end
                for j = 1:obj.ndyn
                    % Qij = { x \in Pi | Ai*x + fi \in Pj }
                    Qij = obj.domain(j).invAffineMap(obj.A{i}, obj.f{i}).intersect(obj.domain(i));
                    if Qij.isFullDim()
                        map.transitions(i, j) = true;
                        map.regions{i, j} = Qij;
                    end
                end
            end
        end

        function [answer, map] = isInvariant(obj, map)
            % Checks whether an autonomous PWA system is invariant
            %
            % Syntax:
            %   answer = pwa.isInvariant()
            %   [answer, map] = pwa.isInvariant()
            %
            % Input:
            %   pwa: autonomous PWA system as an PWASystem object
            %           x^+ = A_i*x+f_i IF x \in R_i
            %
            % Output:
            %   answer: true if the partition of the system is invariant,
            %           false otherwise
            %      map: transition map generated by PWASystem/transitionMap
            
            global MPTOPTIONS
            
            error(obj.rejectArray());
            if obj.nu>0
                error('This method only supports autonomous systems.');
            end
            
            if nargin<2
                % compute the transition map
                fprintf('Computing the transition map...\n');
                map = obj.transitionMap();
                fprintf('...done (%d full-dimensional transitions)\n', nnz(map.transitions));
            end
            
            fprintf('Checking invariance...\n');
            tic;
            answer = true;
            for i = 1:obj.ndyn
                if toc > MPTOPTIONS.report_period
                    fprintf('progress: %d/%d\n', i, obj.ndyn);
                    tic;
                end
                % which parts of Ri enter any other region?
                sources = find(map.transitions(i, :));
                Rsources = [];
                for j = sources
                    Rsources = [Rsources, map.regions{i, j}];
                end
                % are the sources equal to R_i?
                if obj.domain(i)~=Rsources
                    % there exists a subset of R_i which does not enter any
                    % other regions, thus the partition is not invariant
                    answer = false;
                    break
                end
            end
            fprintf('...done\n');
            
        end

        function L = lyapunov(obj, ltype)
            % Constructs a Lyapunov function for an autonomous PWA system
            %
            % Syntax
            %   L = pwa.lyapunov(type)
            %
            % Inputs:
            %   pwa: autonomous PWA system as an PWASystem object
            %           x^+ = A_i*x+f_i IF x \in R_i
            %  type: type of the Lyapunov function ('pwa' or 'pwq')
            %
            % Outputs:
            %     L: PolyUnion object with the Lyapunov function
            
            narginchk(2, 2);
            error(obj.rejectArray());
            if obj.nu>0
                error('This method only supports autonomous systems.');
            end

            switch lower(ltype)
                case 'pwa',
                    lyapfun = @PWALyapFunction;
                case 'pwq'
                    lyapfun = @PWQLyapFunction;
                otherwise
                    error('Unsupported function type "%s".', ltype);
            end
            
            % check invariance and compute the transition map
            [isinv, map] = obj.isInvariant();
            if ~isinv
                error('The system is not invariant.');
            end
            
            % call the corresponding subfunction
            [L, feasible] = feval(lyapfun, obj, map);
            if ~feasible
                error('Infeasible problem.');
            end
        end

    end
    
    methods(Access=private, Hidden)
        % internal APIs
        
        function [L, feasible] = PWALyapFunction(obj, map)
            % PWA Lyapunov function construction for an autonomous PWA system
            %
            % The PWA Lyapunov function is given by
            %   L(x) = L_i'*x + c_i IF x \in R_i
            
            global MPTOPTIONS
            
            % all regions must be bounded
            if ~all(obj.domain.isBounded)
                error('All domains must be bounded.');
            end
            
            % which regions contain the origin?
            origin = zeros(obj.nx, 1);
            containsOrigin = obj.domain.contains(origin);
            % check if the origin is a vertex
            idx = find(containsOrigin);
            for i = 1:length(idx)
                j = idx(i);
                if max(obj.domain(j).A*origin-obj.domain(j).b)<-MPTOPTIONS.abs_tol
                    error('The origin must be a vertex (violated in domain %d).', j);
                end
            end
            
            % prepare variables for each transition
            epsilon = MPTOPTIONS.rel_tol;
            rho = sdpvar(1, 1);
            L = cell(1, obj.ndyn);
            c = cell(1, obj.ndyn);
            for i = 1:obj.ndyn
                L{i} = sdpvar(obj.nx, 1, 'full');
                c{i} = sdpvar(1, 1);
            end
            
            % formulate constraints
            fprintf('Formulating constraints...\n');
            constraints = [ epsilon <= rho <= 1 ];
            tic
            % enforce decay of the Lyapunov function over each transition
            for i = 1:obj.ndyn
                if toc > MPTOPTIONS.report_period
                    fprintf('progress: %d/%d\n', i, obj.ndyn);
                    tic;
                end
                
                % vertices of the i-th region
                % (since we only consider bounded regions, we know the
                % rays are empty)
                Vi = obj.domain(i).V;
                
                % enforce positivity on each vertex
                nx = zeros(size(Vi, 1), 1);
                for j = 1:size(Vi, 1)
                    nx(j) = norm(Vi(j, :)', 1);
                end
                constraints = constraints + ...
                    [ epsilon*nx <= Vi*L{i} + ones(size(nx))*c{i} ];
                
                if containsOrigin(i)
                    % the constant term is zero in all regions which
                    % contain the origin
                    constraints = constraints + [ c{i} == 0 ];
                end
                
                % enforce decrease
                targets = find(map.transitions(i, :));
                for j = targets
                    Vi = map.regions{i, j}.V;
                    nx = zeros(size(Vi, 1), 1);
                    PWAi = [];
                    PWAj = [];
                    for vi = 1:size(Vi, 1)
                        x = Vi(vi, :)';
                        xp = obj.A{i}*x + obj.f{i};
                        PWAi = [PWAi; L{i}'*x + c{i}];
                        PWAj = [PWAj; L{j}'*xp + c{j}];
                        nx(vi) = norm(x, 1);
                    end
                    constraints = constraints + [ rho*nx <= PWAi - PWAj ];
                end
            end
            
            % solve the problem
            fprintf('Solving...\n');
            options = sdpsettings;
            % maximize the decay
            solution = solvesdp(constraints, -rho, options);
            fprintf('...done\n');
            
            % check the solution
            fprintf('\n');
            feasible = true;
            if(solution.problem==4)
                res = checkset(constraints);
                if min(res)>0,
                    fprintf('Numerical problems, but solution is feasible\n');
                else
                    fprintf('Numerical problems, residual: %e (should be positive).\n', min(res));
                end
                
            elseif(solution.problem>0)
                feasible = false;
                res = checkset(constraints);
                fprintf('Infeasible problem, residual: %e (should be positive).\n', min(res));
                L = PolyUnion;
                return
                
            elseif(solution.problem<0)
                error(solution.info);
                
            else
                fprintf('Feasible solution found\n');
            end
            
            % extract the solution
            regions = obj.domain;
            for i = 1:obj.ndyn
                fun = AffFunction(double(L{i})', double(c{i}));
                regions(i).addFunction(fun, 'lyapunov');
            end
            L = PolyUnion(regions);

        end

                
        function [L, feasible] = PWQLyapFunction(obj, map)
            % PWQ Lyapunov function construction for an autonomous PWA system
            %
            % The PWQ Lyapunov function is given by
            %   L(x) = x'*Q_i*x + L_i*x + c_i IF x \in R_i
            
            global MPTOPTIONS
            
            error(obj.rejectArray());
            if obj.nu>0
                error('This method only supports autonomous systems.');
            end

            % all domains must be full-dimensional
            if ~all(obj.domain.isFullDim)
                error('All domains must be full-dimensional.');
            end

            % which regions contain the origin?
            containsOrigin = obj.domain.contains(zeros(obj.nx, 1));
            
            % prepare variables for each transition
            epsilon = MPTOPTIONS.rel_tol;
            rho = sdpvar(1,1);
            Q = cell(1, obj.ndyn);
            L = cell(1, obj.ndyn);
            c = cell(1, obj.ndyn);
            for i = 1:obj.ndyn
                Q{i} = sdpvar(obj.nx, obj.nx, 'symmetric');
                L{i} = sdpvar(obj.nx, 1, 'full');
                c{i} = sdpvar(1, 1);
            end
            
            % formulate constraints
            fprintf('Formulating constraints...\n');
            constraints = [ rho <= -epsilon ];
            tic
            % enforce decay of the Lyapunov function over each transition
            for i = 1:obj.ndyn
                if toc > MPTOPTIONS.report_period
                    fprintf('progress: %d/%d\n', i, obj.ndyn);
                    tic;
                end
                
                % enforce decrease
                targets = find(map.transitions(i, :));
                for j = targets
                    % is this a purely quadratic transition?
                    if containsOrigin(i) && containsOrigin(j)
                        % quadratic transition
                        W = obj.A{i}'*Q{j}*obj.A{i}-Q{i}-rho*eye(obj.nx);
                        constraints = constraints + [ W <= 0 ];
                    else
                        % PWQ transition
                        map.regions{i, j}.minHRep();
                        H = map.regions{i, j}.A;
                        K = map.regions{i, j}.b;
                        m = length(K);
                        N = sdpvar(m, m, 'symmetric');
                        constraints = constraints + [ N(:) >= 0 ];
                        RHS = [-H K]'*N*[-H K];
                        dQ = obj.A{i}'*Q{j}*obj.A{i} - Q{i};
                        dL = 2*obj.A{i}'*Q{j}*obj.f{i} + obj.A{i}'*L{j} - L{i};
                        dc = obj.f{i}'*Q{j}*obj.f{i} + c{j} + obj.f{i}'*L{j} - c{i};
                        W = [dQ-rho*eye(obj.nx), 0.5*dL; 0.5*dL', dc]+RHS;
                        constraints = constraints + [ W <= 0 ];
                    end
                end
                
                % enforce positivity
                P = Q{i}-epsilon*eye(obj.nx);
                if containsOrigin(i)
                    % the function is purely quadratic in a region that
                    % contains the origin
                    constraints = constraints + [ P >= 0 ];
                    constraints = constraints + [ L{i}==0; c{i}==0 ];
                else
                    % the function is positive everywhere else
                    obj.domain(i).minHRep();
                    H = obj.domain(i).A;
                    K = obj.domain(i).b;
                    m = length(K);
                    M = sdpvar(m, m, 'symmetric');
                    constraints = constraints + [ M(:) >= 0 ];
                    RHS = [-H K]'*M*[-H K];
                    W = [P, 0.5*L{i}; 0.5*L{i}', c{i}]-RHS;
                    constraints = constraints + [ W >= 0 ];
                end
            end
            fprintf('...done\n');
            
            % solve the problem
            fprintf('Solving...\n');
            options = sdpsettings;
            solution = solvesdp(constraints, [], options);
            fprintf('...done\n');
            
            % check the solution
            fprintf('\n');
            feasible = true;
            if(solution.problem==4)
                res = checkset(constraints);
                if min(res)>0,
                    fprintf('Numerical problems, but solution is feasible\n');
                else
                    fprintf('Numerical problems, residual: %e (should be positive).\n', min(res));
                end
                
            elseif(solution.problem>0)
                feasible = false;
                res = checkset(constraints);
                fprintf('Infeasible problem, residual: %e (should be positive).\n', min(res));
                L = PolyUnion;
                return
                
            elseif(solution.problem<0)
                error(solution.info);
                
            else
                fprintf('Feasible solution found\n');
            end
            
            % extract the solution
            regions = obj.domain;
            for i = 1:obj.ndyn
                fun = QuadFunction(double(Q{i}), double(L{i})', double(c{i}));
                regions(i).addFunction(fun, 'lyapunov');
            end
            L = PolyUnion(regions);
            
        end
        
    end

end
