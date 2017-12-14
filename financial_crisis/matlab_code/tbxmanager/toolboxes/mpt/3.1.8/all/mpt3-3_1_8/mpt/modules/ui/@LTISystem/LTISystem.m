classdef LTISystem < AbstractSystem
%
%  LTISYSTEM: Represents linear time-invariant systems 
%  ====================================================
%  
%  
%  SYNTAX
%  ------
%     
%      sys = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'T_s', Ts)
%    
%  
%  DESCRIPTION
%  -----------
%     This class represents linear time-invariant systems of the form 
%                            x(t+Ts)  = A x(t) + B u(t)           
%                               y(t)  = C x(t) + D u(t)           
%     where x in R^n_x  is the state vector, u in R^n_u  is the vector of inputs, y
%  in R^n_y  is the vector of outputs, and T_srepresents the sampling time.
%    Each LTI system defines following properties: 
%    
%     - A, B: matrices of the state-update equation (read-only) 
%     - C, D: matrices of the output equation (read-only) 
%     - Ts: sampling time (read-only) 
%     - nx, nu, ny: number of states, inputs and outputs (automatically determined,
%     read-only) 
%     - x: specifications of system's states (see help SystemSignal) 
%     - u: specifications of system's inputs (see help SystemSignal) 
%     - y: specifications of system's outputs (see help SystemSignal) 
%    To define an LTI system, provide the list of system's matrices to the
%  constructor:
%     sys = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'Ts', Ts)All inputs, except
%  of the A matrix, can be omitted. In such a case they are set to empty matrices
%  of corresponding dimension. As a consequence, one can easily define autonomous
%  systems x(t+T_s) = A x(t) + f  by calling sys = LTISystem('A', A, 'f', f, 'Ts',
%  Ts). Similarly, to define an LTI system without outputs, call sys =
%  LTISystem('A', A, 'B', B, 'Ts', Ts). If the sampling time is omitted, it is set
%  to Ts=1.
%    Another option to define an LTI system is to import the dynamics from Control
%  toolbox' discre-time state-space objects:
%     sys = LTISystem(ss(A, B, C, D, Ts))Important to remember is that LTI systems
%  carry with them the value of the state vector. The initial value can be set by
%  the sys.initialize(x0) method (see " help LTISystem/initialize"). Value of the
%  internal state can be retrieved by the sys.getStates() method (see " help
%  LTISystem/getStates"). To update the internal state using the system's
%  state-update equation, use the sys.update(u)function (see " help
%  LTISystem/update").
%  
%  SEE ALSO
%  --------
%     PWASystem
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
        A % Matrix of the state-update equation x(t+Ts) = A*x(t)+B*u(t)+f
        B % Matrix of the state-update equation x(t+Ts) = A*x(t)+B*u(t)+f
		f % Matrix of the state-update equation x(t+Ts) = A*x(t)+B*u(t)+f
        C % Matrix of the output equation y(t) = C*x(t)+D*u(t)+g
        D % Matrix of the output equation y(t) = C*x(t)+D*u(t)+g
		g % Matrix of the output equation y(t) = C*x(t)+D*u(t)+g
	end
	
	methods(Hidden)
		% implementation of abstract methods
		
		% no validation in these functions! it was already performed in
		% AbstractSystem/update() and output()

		function [xn, y] = update_equation(obj, x, u)
			% returns the state update and the output
			
			xn = obj.A*x + obj.B*u + obj.f;
			y = obj.C*x + obj.D*u + obj.g;
		end

		function y = output_equation(obj, x, u)
			% output equation
			
			y = obj.C*x + obj.D*u + obj.g;
		end

		function out = has_feedthrough(obj)
			% feedthrough indication. must return true if the system has
			% direct feedthrough, false otherwise
			
			out = (nnz(obj.D)~=0);
        end
    end

    methods(Access=protected)
        
        function out = display_internal(obj)
            % Returns a string information about the system
            
            plural = @(s, n) [num2str(n) ' ' s repmat('s', 1, double(n~=1))];
            out = sprintf('%s with %s, %s, %s', class(obj), ...
					plural('state', obj.nx), ...
					plural('input', obj.nu), ...
					plural('output', obj.ny));
        end
    end
    
    methods
        
        function obj = LTISystem(varargin)
            % Constructor for LTI systems
            %
            % To create an LTI system
            %   x^+ = Ax + Bu + f
            %     y = Cx + Du + g
            % call:
            %   s = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g)
            %
            % To import LTI system defined by a state space object call
            %   s = LTISystem(object)
			%
            % To import from a sysStruct structure:
            %   s = LTISystem(sysStruct)
            
            if nargin == 0
                return
			end
            
			if nargin==1 && isstruct(varargin{1})
				% import from sysStruct
                S = mpt_verifySysStruct(varargin{1});
				if iscell(S.A)
					error('Use PWASystem for PWA systems.');
				end
				if ~isfield(S, 'f')
					S.f = [];
				end
				if ~isfield(S, 'g')
					S.g = [];
				end
				
			elseif nargin==1 && isa(varargin{1}, 'ss') && ...
					isdt(varargin{1})
                % import from state space object
                ss = varargin{1};
				S.A = ss.A;
				S.B = ss.B;
				S.C = ss.C;
				S.D = ss.D;
				S.f = [];
				S.g = [];
				S.Ts = ss.Ts;
				
			elseif nargin==1
				
				error('Can only import from discrete-time state-space objects.');
			
			else
				
				ip = inputParser;
				ip.KeepUnmatched = false;
				ip.addParamValue('A', [], @validate_realmatrix);
				ip.addParamValue('B', [], @validate_realmatrix);
				ip.addParamValue('f', [], @validate_realmatrix);
				ip.addParamValue('C', [], @validate_realmatrix);
				ip.addParamValue('D', [], @validate_realmatrix);
				ip.addParamValue('g', [], @validate_realmatrix);
				ip.addParamValue('domain', [], @validate_polyhedron);
				ip.addParamValue('Ts', 1, @(x) isa(x, 'double') && numel(x)==1 && x>=0);
				ip.parse(varargin{:});
				S = ip.Results;
			end
			if ~isfield(S, 'domain')
				S.domain = [];
			end
                
            obj.A = S.A;
            obj.B = S.B;
            obj.C = S.C;
            obj.D = S.D;
            try
                obj.Ts = S.Ts;
			end
            
            obj.nx = size(obj.A, 2);
            obj.nu = size(obj.B, 2);
            obj.ny = size(obj.C, 1);

			if isempty(obj.B)
				obj.B = zeros(obj.nx, obj.nu);
			end
			if isempty(obj.C)
				obj.C = zeros(obj.ny, obj.nx);
			end
			if isempty(obj.D)
				obj.D = zeros(obj.ny, obj.nu);
			end
			if isempty(S.f)
				S.f = zeros(obj.nx, 1);
			end
			if isempty(S.g)
				S.g = zeros(obj.ny, 1);
			end
			obj.f = S.f;
			obj.g = S.g;

			if size(obj.A, 1)~=size(obj.A, 2)
				error('The "A" matrix must be a %dx%d matrix.', obj.nx, obj.nx);
			end
			if size(obj.A, 1)~=size(obj.B, 1)
				error('The "B" matrix must have %d rows.', obj.nx);
            end
            if size(obj.C, 2)~=size(obj.A, 2)
                error('The "C" matrix must have %d columns.', obj.nx);
            end
			if size(obj.C, 1)~=size(obj.D, 1)
				error('The "D" matrix must have %d rows.', obj.ny);
			end
			if size(obj.B, 2)~=size(obj.D, 2)
				error('The "B" and "D" matrices must have the same number of columns.');
			end
			if size(obj.f, 1)~=obj.nx || size(obj.f, 2)~=1
				error('"f" must be a %dx1 vector.', obj.nx);
			end
			if size(obj.g, 1)~=obj.ny || size(obj.g, 2)~=1
				error('"g" must be a %dx1 vector.', obj.ny);
			end

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
            % Converts the LTI system into a state-space object
            
			if nnz(obj.f)>0 || nnz(obj.g)>0
				error('System has affine terms, cannot convert to SS.');
			end
			if isempty(obj.Ts)
				error('Sampling time must be specified.');
			end
            out = ss(obj.A, obj.B, obj.C, obj.D, obj.Ts);
        end
        
		function obj = setDomain(obj, type, domain)
			% Sets domain of the system
			
			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			narginchk(3, 3);
			if ~ischar(type)
				error('The first input must be a string.');
			end
			if ~isa(domain, 'Polyhedron')
				error('The second input must be a polyhedron.');
			end

			function check_dimension(domain, required_dimension)
				if domain.Dim ~= required_dimension
					error('The domain must have dimension %d.', required_dimension);
				end
			end

			% this anonymous function constructs an infinity box of
			% dimension 'n', which emulates R^n
			infbox = @(n) Polyhedron('lb', -MPTOPTIONS.infbound*ones(n, 1), 'ub', MPTOPTIONS.infbound*ones(n, 1));

			domain.minHRep();
			switch type
				case 'x',
					% check and extend the domain to x-u space
					check_dimension(domain, obj.nx);
					obj.domain = domain*infbox(obj.nu);
					obj.domainx = domain;

				case 'u'
					% check and extend the domain to x-u space
					check_dimension(domain, obj.nu);
					obj.domain = infbox(obj.nx)*domain;
					
					% TODO: remove the hard-coded 'fourier' method once
					% Polyhedron/projection is reliable
					obj.domainx = obj.domain.projection(1:obj.nx, 'fourier');
					
				case 'xu',
					% x-u domain provided, just check
					check_dimension(domain, obj.nx+obj.nu);
					obj.domain = domain;
					
					% TODO: remove the hard-coded 'fourier' method once
					% Polyhedron/projection is reliable
					obj.domainx = obj.domain.projection(1:obj.nx, 'fourier');
					
				otherwise
					error('Unrecognized option "%s".', type);
			end
		end
        
		function K = LQRGain(obj)
			% Returns the LQR gain u=K*x
			
			if nnz(obj.f)>0 || nnz(obj.g)>0
				error('This function does not support affine systems.');
			end
			error(obj.assert_has_xu_penalties);
			K = -dlqr(obj.A, obj.B, obj.x.penalty.weight, obj.u.penalty.weight);
		end
		
		function P = LQRPenalty(obj)
			% Returns the LQR penalty

			if nnz(obj.f)>0 || nnz(obj.g)>0
				error('This function does not support affine systems.');
			end
			error(obj.assert_has_xu_penalties);
			[~, Q] = dlqr(obj.A, obj.B,	obj.x.penalty.weight, obj.u.penalty.weight);
			P = QuadFunction(Q);
		end
		
		function S = LQRSet(obj)
			% Returns the LQR invariant set

			if nnz(obj.f)>0 || nnz(obj.g)>0
				error('This function does not support affine systems.');
			end
			error(obj.assert_has_xu_penalties);
			CL = ClosedLoop(LQRController(obj), obj);
			S = CL.toSystem().invariantSet();
		end
		
		function C = stabilizingController(obj, type)
			% Computes a stabilizing controller

			if nnz(obj.f)>0 || nnz(obj.g)>0
				error('This function does not support affine systems.');
			end

			if nargin<2
				type = 'pwq';
			end
			
			if isequal(type, 'pwq')
				L = LQRController(obj).toInvariant();
				C = EMPCController(L);
					
			else
				error('not yet implemented.');
			end			
		end
		
		function [S, SN] = reachableSet(obj, varargin)
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
			
			backwards = lower(options.direction(1))=='b';
			X = options.X;
			U = options.U;
			SN = {}; % reachable sets at the i-th step

			% lift the system to x-u space
			Acl = [obj.A, obj.B; zeros(obj.nu, obj.nx) eye(obj.nu)];
			fcl = [obj.f; zeros(obj.nu, 1)];
			
			% compute the N-step reachable set
			for i = 1:options.N
				SN{i} = [];
				for j = 1:numel(X)
					% we support multiple state domains
					%
					% TODO: replace the for-loop by Polyhedron/forEach

					if X(j).isEmptySet()
						continue
					end

					D = X(j)*U;
					if backwards
						XUn = D.invAffineMap(Acl, fcl);
					else
						XUn = Acl*D+fcl;
					end
					if backwards
						XUn = XUn.intersect(obj.domain);
					end

					if XUn.isEmptySet()
						continue
                    end
                    if obj.nu > 0
                        Xn = XUn.projection(1:obj.nx);
                    else
                        Xn = XUn;
                    end
					Xn.minHRep();
					SN{i} = [SN{i}, Xn];
				end
				X = SN{i};
			end
			S = SN{end};
            if isempty(S)
                S = Polyhedron.emptySet(obj.nx);
            end
		end

		function [X, dynamics] = invariantSet(obj, varargin)
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
            ip.addParamValue('U', ...
				obj.u.boundsToPolyhedron(), ...
				@validate_polyhedron);
            ip.addParamValue('X', ...
				obj.x.boundsToPolyhedron(), ...
				@validate_polyhedron);
			ip.parse(varargin{:});
			options = ip.Results;
			if obj.nu>0 && options.U.isEmptySet()
				error('Input constraints must not be empty.');
			end
			if obj.nu>0 && options.U.Dim~=obj.nu
				error('Input constraints must be a polyhedron in %dD.', obj.nu);
			end
			if options.X.Dim~=obj.nx
				error('State constraints must be a polyhedron in %dD.', obj.nx);
			end

			Xo = options.X;
			U = options.U;
			converged = false;
			for i = 1:options.maxIterations
				fprintf('Iteration %d...\n', i);
				X = obj.reachableSet('X', Xo, 'U', U, ...
					'direction', 'backward');
				X = X.intersect(Xo).minHRep();
				if X==Xo
					converged = true;
					break
				else
					Xo = X;
				end
			end
			if ~converged
				warning('Computation finished without convergence.');
			end
			dynamics = ones(1, length(X));
		end
                
		function C = constraints(obj)
			% Convert LTI model into YALMIP constraints
			
			% constraints on variables
			C = constraints@AbstractSystem(obj);
			
			% add the LTI dynamics constraints
			x = obj.x.var;
			u = obj.u.var;
			y = obj.y.var;
			for k = 1:obj.Internal.system.N
				if obj.nx > 0
					C = C + [ x(:, k+1) == obj.A*x(:, k) + obj.B*u(:, k) + obj.f];
				end
				if obj.ny > 0
					C = C + [ y(:, k) == obj.C*x(:, k) + obj.D*u(:, k) + obj.g];
				end
			end
		end
	end

end
