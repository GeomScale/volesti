classdef SystemSignal < FilterBehavior & IterableBehavior
%
%  SYSTEMSIGNAL: Class representing variables of a dynamical system 
%  =================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      signal = SystemSignal(n)
%    
%  
%  DESCRIPTION
%  -----------
%     The SystemSignal class represents prediction of variables (e.g., system
%  states, inputs, or outputs) of a dynamical system in an MPC setup.
%    Signals are created by signal = SystemSignal(n), where n is the dimension of
%  the variable. Once constructed, signals can be assigned names by signal.name =
%  'myname'. Each signal has following properties available for read/write access
%  by default: 
%    
%     - signal.min: lower bound on the variable in MPC problems 
%     - signal.max: upper bound on the variable in MPC problems 
%     - signal.penalty: penalization of the variable in the MPC cost function given
%     as Function object 
%    Many additional properties can be set by the concept of filters. In short,
%  filters are dynamical properties which can be added to a signal on demand. A
%  filter is added by calling signal.with('filter_name'). List of available filters
%  can be obtained by calling methods SystemSignal and looking for methods prefixed
%  by the filter_ string. Filters can be removed by calling
%  signal.without('filter_name'). To list filters added to a particular signal, use
%  the signal.listFilters() method.
%  
%  INPUT
%  -----
%     
%        
%          n Dimension of the variable                
%            Class: double                            
%              
%  
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
 
 
    properties
        name % Name of the signal
    end
    properties(Hidden)
		type % Type of the signal (input, output, state)
		init % Initial value of the signal
        userData % User-defined value
        previous_value = [] % Previous value of the signal
    end
    properties(SetAccess = private, Hidden)
        var % YALMIP variable representing the signal
        n = 0 % Dimension of the signal
        N % Length of the signal
	end
	
    methods
        
        function obj = SystemSignal(n)
            % constructor:
            %   s = SystemSignal(n)
            %           n: number of variables

            if nargin == 0
                return
			end
            
            obj.n = n;
            obj.N = 1;

			% each signal has the min/max properties by default
			obj.with('min');
			obj.with('max');
			obj.with('penalty');
			
		end

		function s = saveobj(obj)
			% save method

			% we need to work on a copy of the object, since we are going
			% to remove from it filters and sdpvar objects
			s = copy(obj);
			s.saveSdpvarValue();
			% remove the sdpvar
			s.uninstantiate();
			s.saveAllFilters();
		end
		
		function out = join(obj)
			% Joins several signals into a single one
			
			M = numel(obj);
			if M<2
				out = obj;
				return
			end
			
			error('Not yet implemented.');
			
			n = zeros(1, M); instantiated = zeros(1, M); N = zeros(1, M);
			for i = 1:M
				n(i) = obj(i).n;
				instantiated(i) = ~isempty(obj(i).var);
				N(i) = obj(i).N;
			end
			out = SystemSignal(sum(n));
			
			if all(instantiated)
				if all(diff(N)==0)
					% instantiate
					out.instantiate(N(1));
				else
					error('Cannot join signals with different prediction horizons.');
				end
			end
			
		end
		
        function out = isKind(obj, kind)
            % Returns the kind of the variable ('x', 'u', 'y', ...)
            
            out = isfield(obj.userData, 'kind') && ...
                ((ischar(kind) && isequal(obj.userData.kind, kind)) || ...
                (iscell(kind) && ismember(obj.userData.kind, kind)));
		end
		
		function setKind(obj, kind)
			% Sets the kind of a signal
			
			obj.userData.kind = kind;
		end
        
        function obj = initialize(obj, init)
            % Sets the initial value of the signal
            
            obj.init = init;
		end
        
        
		function obj = uninstantiate(obj)
			% removes YALMIP's representation of the signal
			
			if isa(obj.var, 'sdpvar')
				obj.var = [];
				% notify filters that they should remove any
				% self-introduced YALMIP variables
				obj.applyFilters('uninstantiate');
			end
		end
		
        function obj = instantiate(obj, N)
            % instantiate variable on a prediction horizon
            
            if obj.isKind('x')
                M = N + 1;  % we need one more state variable
            else
                M = N;
            end
            obj.N = M;

			% create YALMIP variables
			obj.var = sdpvar(obj.n, M, 'full');

			obj.applyFilters('instantiate');

            %             N_required = max(0, M-size(obj.var, 2));
            %             if N_required > 0
            %                 if obj.isbin
            %                     obj.var = [obj.var binvar(obj.n, N_required, 'full')];
            %                 else
            %                     obj.var = [obj.var sdpvar(obj.n, N_required, 'full')];
            %                 end
            %             end
            
        end

        function C = constraints(obj)
            % Convert variable into YALMIP constraints

            % add custom filters
            C = obj.applyFilters('constraints');
            
        end
        
        function out = objective(obj)
            % Convert variable into objective

            % add custom filters, e.g. delta penalties
            out = obj.applyFilters('objective');
            
        end
        
        function out = value(obj, k)
            % Returns optimized value of the signal at time step "k"
            
            out = double(obj.var(:, 1:obj.N));
            if nargin == 2 && ~isempty(out);
                out = out(:, k);
            end
        end
        
        function plot(obj)
            % Plots optimal values of the signal
            
			if isempty(obj.var)
				error('Cannot plot signal which was not optimized.');
			end
			v = obj.value();
			if any(isnan(v))
				error('Cannot plot signal which was not optimized.');
			end
            stairs(0:obj.N, [v obj.value(obj.N)]', 'LineWidth', 2);
            grid on
            
            % enlarge axis to improve readability
            min_v = min(min(v, [], 2));
            max_v = max(max(v, [], 2));
            d = max_v - min_v;
            min_v = min_v - 0.1*d;
            max_v = max_v + 0.1*d;
            axis([0 obj.N min_v-1e-10 max_v+1e-10]);
		end
		
		function out = isConstrained(obj)
			% Returns true if signal is subject to constraints
			
			out = any(~isinf(obj.max)) || any(~isinf(obj.min)) || ...
				obj.hasFilter({'setConstraint', 'terminalSet'});
		end
		
		function out = boundsToPolyhedron(obj)
			% Converts signal constraints to a polyhedron object
			
			if obj.n==0
				out = Polyhedron;
			else
				out = Polyhedron([eye(obj.n); -eye(obj.n)], ...
					[sanitize_inf(obj.max); -sanitize_inf(obj.min)]);
                if obj.hasFilter('setConstraint')
                    out = out & obj.setConstraint;
                end
			end
		end
	end
	
	methods(Access = protected)
		
% 		function new = copyElement(obj)
% 			% Copy constructor
% 			
% 			new = copyElement@FilterBehavior(obj);
% 			
% 			if ~isempty(obj.var)
% 				% re-instantiate variables
% 				
% 				if obj.isKind('x')
% 					% accommodate for the fact that
% 					% SystemSignal/instantiate automatically adds one more
% 					% step for state variables
% 					new.instantiate(obj.N-1);
% 				else
% 					new.instantiate(obj.N);
% 				end
% 				
% 				% re-assign values
% 				%assign(new.var, double(obj.var));
% 			end
% 			
% 		end

	end

	methods(Hidden)

		function msg = validatePenalty(obj, P)
			% validates penalty P

			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			if ~isa(P, 'Function')
				msg = 'Input must be a Function object.';
				return
			else
				msg = '';
			end
			
			if isempty(P.findprop('weight'))
				% "P" is a Function object, no validation is possible, use
				% at your own risk
				return
			end
			
			Q = P.weight;
			check_definiteness = isa(P, 'QuadFunction');
				
			if size(Q, 2) ~= obj.n
				msg = sprintf('The weighting matrix must have %d column(s).', obj.n);
			elseif check_definiteness
				if size(Q, 1) ~= size(Q, 2)
					msg = 'The weighting matrix must be square.';
				elseif obj.isKind('u') && min(eig(Q)) < MPTOPTIONS.abs_tol
					msg = 'The weighting matrix must be positive definite.';
				elseif min(eig(Q)) < -MPTOPTIONS.abs_tol
					msg = 'The weighting matrix must be positive semi-definite.';
				end
			end
		end
		
		function obj = saveSdpvarValue(obj)
			% saves the current value of the sdpvar to
			% obj.Internal.save.value

			if isa(obj.var, 'sdpvar')
				obj.Internal.save.value = double(obj.var);
			end
		end

		function obj = loadSdpvarValue(obj)
			% loads value of an sdpvar stored in
			% obj.Internal.save.value

			if isfield(obj.Internal, 'save') && ...
					isfield(obj.Internal.save, 'value')
				value = obj.Internal.save.value;
				if ~isa(obj.var, 'sdpvar')
					% re-instantiate first
					obj.var = sdpvar(size(value, 1), size(value, 2), 'full');
				end
				assign(obj.var, value);
			end
		end

		function initializeFromRecord(obj, data)
            
            if obj.isKind('x')
                obj.init = data(:, end);
                if ~isempty(obj.var)
                    assign(obj.var(1), data(:, end));
                end
            end
            if ~isempty(data)
                obj.previous_value = data(:, end);
            end
        end

        function data = record(obj, data)
            
            if nargin == 1
                if obj.isKind('x')
                    data = obj.init;
                else
                    data = [];
                end
            else
                if obj.isKind('x')
                    data = [data obj.value(2)];
                else
                    data = [data obj.value(1)];
                end
            end
		end

	end
end
            
