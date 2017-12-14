classdef MPCController < AbstractController
%
%  MPCCONTROLLER: MPC controller with on-line optimization 
%  ========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ctrl = MPCController(model)
%      ctrl = MPCController(model, horizon)
%    
%  
%  DESCRIPTION
%  -----------
%     ctrl = MPCController(model) constructs an MPC optimization problem by using
%  the object model as the prediction model.
%    The basic optimization problem takes the following form: 
%                                        N-1                         
%                                        --                          
%                                        \                           
%                           min  J(x ) + /   J(x , u , y )           
%                            U      N    --     k   k   k            
%                                        k=0                         
%                                                                    
%                          s.t.  x    = f(x , u )                    
%                                 k+1      k   k                     
%                                                                    
%                                y  = g(x , u )                      
%                                 k      k   k                       
%                                                                    
%                                x    <= x  <= x                     
%                                 min     k     max                  
%                                                                    
%                                u    <= u  <= u                     
%                                 min     k     max                  
%                                                                    
%                                y    <= y  <= y                     
%                                 min     k     max                  
%     where N  is the prediction horizon, x_k, u_k, y_k  are the state, input, and
%  output predictions, respectively, J(x_N)  denotes the terminal cost, J(x_k, u_k,
%  y_k)  is the stage cost, f(x_k, u_k)represents the state-update equation, and
%  g(x_k, u_k)  is the output equation.
%    The terminal cost is given by J(x_N) = || Q x_N ||_p, where Q  is the
%  weighting matrix and p  is the type of the norm (if p=2, then J(x_N) = x_N^T Q
%  x_N, see " help Function" for more information). The type of the norm as well as
%  the weighting matrix are specified by the terminalPenalty property of the state
%  signal (see " help SystemSignal/filter_terminalPenalty").
%    The stage cost is represented by J(x_k, u_k, y_k) = || Q_x x_k ||_p_x + || Q_u
%  u_k ||_p_u + || Q_y y_k ||_p_y, where Q_i  are weighting matrices and p_i 
%  represent types of the norms. Properties of the state penalization are taken
%  from model.x.penalty (see " help SystemSignal/filter_penalty"). Similarly,
%  properties of input and output penalizations are extracted from model.u.penalty
%  and model.y.penalty, respectively. Note that different norms/weights can be
%  assigned to various signals. Hence you can have a quadratic penalization of the
%  state, but use a 1-norm penalty on the inputs.
%    If a certain signal (state, input, and/or output) has the reference filter
%  enabled (see " help SystemSignal/filter_reference"), the terminal and stage
%  costs become J(x_N) = || Q (x_N - x_ref) ||_p  and J(x_k, u_k, y_k) = || Q_x
%  (x_k - x_ref ||_p_x + || Q_u (u_k - u_ref) ||_p_u + || Q_y (y_k - y_ref) ||_p_y.
%  If some signal does not have the reference property enabled, we assume a zero
%  reference instead.
%    By default, min/max bounds on state, input, and output variables are added to
%  the MPC setup. The bounds stored in model.x.min, model.x.max, model.u.min,
%  model.u.max, model.y.min, and model.y.max are used. Note that other types of
%  constraints can be added. These include polyhedral constraints (see " help
%  SystemSignal/filter_setConstraint"), soft constraints (see " help
%  SystemSignal/filter_softMax" and " help SystemSignal/filter_softMin"), or move
%  blocking (see " help SystemSignal/filter_block").
%    The basic way to create an MPC controller is to call
%     controller = MPCController(system)
%    where system represents the dynamical system to be used for predictions (can
%  be an instance of the LTISystem, PWASystem, or MLDSystem classes). Then you can
%  further fine-tune the prediction model by adding constraints and/or modifying
%  penalties. This can be done by accessing the modelproperty of the controller
%  object:
%     controller.model.x.min = [-5; -5]
%    You can also specify the prediction horizon directly from the command line:
%     controller.N = 3
%    Once the controller is fully specified, you can use the evaluate method to
%  obtain the optimal sequence of inputs, (see " help MPCController/evaluate"), or
%  compute the explicit form of the controller by calling the toExplicitmethod (see
%  " help MPCController/toExplicit").
%  
%  INPUT
%  -----
%     
%        
%          model Any MPT3 system ( LTISystem, PWASystem,  
%                MLDSystem)                               
%                Class: AbstractSystem                    
%          N     Prediction horizon                       
%                Class: Double                            
%                  
%  
%  
%  OUTPUT
%  ------
%     
%        
%          controller Instance of the MPCController class.     
%                       
%  
%  
%  SEE ALSO
%  --------
%     EMPCController
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
 
 
    methods
		
		function out = getName(obj)
			out = 'MPC controller';
		end
		
		function out = isExplicit(obj)
			out = false;
		end

        function obj = MPCController(varargin)
            % Constructor:
            %   ctrl = Controller(model, N)
			%         model: prediction model
			%             N: prediction horizon

			%obj.addlistener('model', 'PreGet', @obj.modelPreGetEvent);
			%obj.addlistener('model', 'PostGet', @obj.modelPostGetEvent);
			%obj.addlistener('model', 'PreSet', @obj.modelPreSetEvent); 
			%obj.addlistener('model', 'PostSet', @obj.modelPostSetEvent); 
			obj.addlistener('N', 'PostSet', @obj.NPostSetEvent);

			if nargin==0
				return
			end

			obj.importUserData(varargin{:});
			
			% Let's not construct the optimizer here, since it is expected
			% that users will modify the controller often.
			% 			if ~isempty(obj.N)
			% 				% If prediction horizon was provided, construct the
			% 				% optimizer immediately. Otherwise it will be constructed
			% 				% once evaluate() is called for the first time (also
			% 				% applies if the prediction horizon is changed afterwards).
			% 				obj.construct();
			% 			end
		end
		
        function EMPC = toExplicit(obj)
            % Computes the explicit controller
            
			% make sure we have the prediction horizon available
			error(obj.assert_controllerparams_defined);
			
            EMPC = EMPCController;
			EMPC.model = obj.model;
			EMPC.N = obj.N;
			
			% propagate custom YALMIP setup if provided
			if ~isempty(obj.yalmipData)
				% use YALMIP data to construct the explicit solution
				EMPC.fromYALMIP(obj.yalmipData);
			else
				% construct the explicit solution
				EMPC.construct();
			end
		end
        		
		function obj = construct(obj, sdpopt)
			% Converts the MPC problem into a YALMIP's optimizer object
            %
            % Syntax:
            %   ctrl.construct()
            %   ctrl.construct(sdpopt) -- uses the specified sdpsettings

			% make sure we have the prediction horizon available
			error(obj.assert_controllerparams_defined);
			
			if isempty(obj.yalmipData)
				Y = obj.toYALMIP();
			else
				Y = obj.yalmipData;
            end
			
            if nargin==2
                Y.internal.sdpsettings = sdpopt;
            end
			% construct YALMIP's optimizer object with the first state as
			% the initial condition
			obj.optimizer = optimizer(Y.constraints, Y.objective, ...
				Y.internal.sdpsettings, ...
				Y.internal.parameters, Y.internal.requested);
			
			obj.markAsUnmodified();
		end
		
        function [u, feasible, openloop] = evaluate(obj, xinit, varargin)
            % Solves the MPC optimization problem for a given initial
            % condition
			%
			% u = controller.evaluate(x0) solves the MPC optimization
			% problem using x0 as the initial condition and returns the
			% first element of the optimal sequence of control inputs. If
			% the problem is infeasible, "u" will be NaN.
			%
			% [u, feasible] = controller.evaluate(x0) also returns a
			% boolean flag indicating whether the optimization problem had
			% a feasible solution.
			%
			% [u, feasible, openloop] = controller.evaluate(x0) also
			% returns the open-loop predictions of states, inputs and
			% outputs in "openloop.X", "openloop.Y", and "openloop.Y",
			% respectively. Value of the optimal cost is returned in
			% openloop.cost. 
			%
			% u = controller.evaluate(x0, 'x.reference', ref, 'u.prev', u0)
			% includes "ref" and "u0" into the vector of initial
			% conditions.

			% make sure we have the prediction horizon available
			error(obj.assert_controllerparams_defined);
			
			error(validate_vector(xinit, obj.nx, 'initial state'));
			
			if obj.wasModified() || isempty(obj.optimizer)
				% Re-generate the optimizer if the problem setup was
				% changed
				obj.construct();
			end

			% assemble the vector of initial conditions. Include any
			% variables that were declared by filters as those which need
			% proper initialization.
			xinit = obj.parse_xinit(xinit, varargin{:});
			
			% use a pre-constructed optimizer object for faster
			% evaluation
			[U, status] = obj.optimizer{xinit};
            % these statuses indicate a feasible solution 
            % (from "help yalmiperror"):
            %  0: Successfully solved
            %  3: Maximum #iterations or time-limit exceeded
            %  4: Numerical problems
            %  5: Lack of progress
            feasible = ismember(status, [0, 3, 4, 5]);
			if ~feasible
				J = Inf;
				u = NaN(obj.nu, 1);
				U = NaN(size(U));
			elseif isempty(U)
				error('Ooops, something went wrong. Please report to mpt@control.ee.ethz.ch');
			else
				J = U(1);
				u = U(2:obj.nu+1);
			end

			if nargout==3
				% also return the open-loop profiles
				%
				% TODO: also return the 'd' and 'z' variables
				openloop.cost = J;
				U = U(2:end);
				openloop.U = reshape(U(1:obj.nu*obj.N), [obj.nu obj.N]);
				U = U(obj.nu*obj.N+1:end);
				openloop.X = reshape(U(1:obj.nx*(obj.N+1)), [obj.nx, obj.N+1]);
				U = U(obj.nx*(obj.N+1)+1:end);
				openloop.Y = reshape(U(1:obj.model.ny*obj.N), [obj.model.ny, obj.N]);
			end
		end
		
		function new = saveobj(obj)
			% save method

			new = saveobj@AbstractController(obj);

			% YALMIP's @optimizer objects cannot be saved
			new.optimizer = [];
			
        end
        
        function Matrices = getMatrices(obj)
            % Returns matrices of the parametric problem
            %
            %   Matrices = ctrl.getMatrices()
            %
            % returns the matrices of the pQP/pLP formulation of the MPC
            % optimization problem. If Matrices.qp=1, the problem is
            % formulated as
            %   min_U 1/2 U'*H*U + x'*F*U + Cf*U + x'*Y*x + Cx*x + Cc
            %    s.t. G*U <= W + E*x
            %
            % If the problem is a pLP (identified by Matrices.qp=0), then:
            %   min_U H*U + Cx*x + Cc
            %    s.t. G*U <= W + E*x
            %
            % Here, U is the open-loop sequence of optimal control inputs
            % and "x" is the vector of parameters.
            
            % make sure we have the prediction horizon available
			error(obj.assert_controllerparams_defined);
            
            % export the problem to YALMIP
            Ydata = obj.toYALMIP();
            
            Matrices = mpt_yalmip2mpt(Ydata.constraints, Ydata.objective, ...
                Ydata.internal.parameters(:), Ydata.variables.u(:));
        end
        
        function f = toFiordos(obj)
            % Exports the MPC controller to FiOrdOs
            %
            %   ctrl.toFiordos()
            
            Y = obj.toYALMIP();
            opt = optimizer(Y.constraints, Y.objective, [], ...
                Y.internal.parameters, Y.variables.u(:));
            f = fiordos(opt);
        end
	end

	methods(Access=protected)
		
		function status = wasModified(obj)
			% Returns true if either the prediction horizon or the problem
			% properties (e.g. constraints) were changed by the user
			
			status = obj.wasModified@AbstractController || ...
				any(cellfun(@any, obj.model.forEachComponent('all', 'wasModified')));
		end
		
		function obj = markAsUnmodified(obj)
			% Marks the object as unmodified
			
			obj.markAsUnmodified@AbstractController;
			obj.model.forEachComponent('all', 'markAsUnmodified');
		end

	end
	
	methods(Access=protected, Hidden)
		
		
		function modelPostSetEvent(obj, source, event)
		end

		function modelPreSetEvent(obj, source, event)
		end

		function modelPreGetEvent(obj, source, event)
		end

		function modelPostGetEvent(obj, source, event)
		end
		
		function NPreSetEvent(obj, source, event)
		end

		function NPostSetEvent(obj, source, event)
			% triggered after the prediction horizon was changed
			
			obj.markAsModified();
		end

	end

	methods(Static, Hidden)
		
		function checkOptimizer(optimizer)
			% validate that the optimizer is either [] or an instance of
			% the @optimizer class
			
			if ~((isa(optimizer, 'double') && isempty(optimizer)) || ...
					isa(optimizer, 'optimizer'))
				error('The optimizer must be an instance of the @optimizer class.');
			end
		end
	end
end
