classdef SFController < EMPCController
    % Class representing state-feedback controllers
    %
    % Constructor:
	%   ctrl = SFController(model, gain)
    
    methods
		
		function out = getName(obj)
			out = 'State feedback controller';
		end

		function obj = construct(obj, K, P_N)
            % Constructs the explicit solution

            if isempty(P_N)
                Q = QuadFunction(zeros(obj.model.nx));
            else
                Q = QuadFunction(P_N);
            end
			P = obj.model.domainx;
			P.addFunction(AffFunction(K, zeros(obj.model.nu, 1)), 'primal');
			P.addFunction(QuadFunction(Q.weight), 'obj');
			obj.optimizer = PolyUnion('Set', P, 'Bounded', true, ...
				'FullDim', true, 'Convex', true);
		end
		
		function Y = toYALMIP(obj)
			error('This controller cannot be exported to YALMIP.');
		end
		
		function obj = fromYALMIP(obj, Y)
			error('This controller cannot be imported to YALMIP.');
		end

        function obj = SFController(model, K, P)
            % Constructor:
			%
			%   ctrl = SFController(model)
            
            if nargin==0
                return
            end
            if nargin<3
                P = [];
            end
            assert(isa(model, 'AbstractSystem'), 'The first input must be a system model.');
            assert(isa(K, 'double'), 'The gain must be a double matrix.');
			obj.importUserData(model);
            if model.nu~=size(K, 1)
                error('The gain must have %d row(s).', model.nu);
            end
            if model.nx~=size(K, 2)
                error('The gain must have %d column(s).', model.nx);
            end
			obj.N = 1;
			if ~isobject(obj.optimizer)
				% construct the explicit solution
				obj.construct(K, P);
			end
		end
		

    end
end
