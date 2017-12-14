classdef mptctrl < AbstractController
    % Compatibility class representing MPT2 controllers

	methods
		function obj = mptctrl()
		end
	end

	methods(Static)
		function out = loadobj(ctrl)
			% load method
			
			% we simply convert explicit MPC to EMPCController objects and
			% on-line MPC to MPCController objects.
			
			fprintf('Loading an MPT2 controller...\n');
			[sysStruct, probStruct] = mpt_verifySysProb(ctrl.details.origSysStruct, ...
				ctrl.details.origProbStruct);
			if isinf(probStruct.N)
				% TODO: revisit once we start supporting CITOC problems
				fprintf('Warning: Infinite prediction horizon changed to N=1.\n');
				probStruct.N = 1;
			end
			mintime = (probStruct.subopt_lev==1);
			probStruct.subopt_lev=0;

			model = mpt_import(sysStruct, probStruct);
			if isequal(ctrl.type, 'explicit')
				% TODO: detect convexity of the partition (not easy, since
				% even for LTI systems we could get non-convex partitions
				% e.g. when binary inputs were used)
				ctrl.convex = 0;
				if mintime
					out = EMinTimeController(mpt_mpsol2pu(ctrl));
				else
					out = EMPCController(mpt_mpsol2pu(ctrl));
				end
			else
				out = MPCController;
			end
			out.model = model;
			out.N = probStruct.N;
		end
	end
end
