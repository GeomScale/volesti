classdef LQRController < EMPCController
    % Class representing LQR state feedback
    %
    % Constructor:
	%   ctrl = LQRController(model)
    
    methods
		
		function out = getName(obj)
			out = 'LQR controller';
		end

		function obj = construct(obj)
            % Constructs the explicit solution

			K = obj.model.LQRGain();
			Q = obj.model.LQRPenalty();
			
			P = obj.model.domainx;
			P.addFunction(AffFunction(K, zeros(obj.model.nu, 1)), 'primal');
			P.addFunction(QuadFunction(Q.weight), 'obj');
			
			obj.optimizer = PolyUnion('Set', P, 'Bounded', true, ...
				'FullDim', true, 'Convex', true);
		end
		
		function Y = toYALMIP(obj)
			error('LQR setups cannot be exported to YALMIP.');
		end
		
		function obj = fromYALMIP(obj, Y)
			error('LQR setups cannot be imported to YALMIP.');
		end

        function obj = LQRController(varargin)
            % Constructor:
			%
			%   ctrl = LQRController(model)
            
			if nargin==0
				return
			end
			
			obj.importUserData(varargin{:});
			obj.N = 1;
			if ~isobject(obj.optimizer)
				% construct the explicit solution
				obj.construct();
			end
		end
		

    end
end
