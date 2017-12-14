function ctrl = mpt_control(sysStruct, probStruct, type, varargin)
% To import MPC properties from sysStruct and probStruct:
%   model = mpt_import(sysStruct, probStruct)
%
% To design an on-line MPC controller
%   online_ctrl = MPCController(model, prediction_horizon)
%
% To compute the explicit representation of the controller:
%   exp_ctrl = online_ctrl.toExplicit()
%
% For more information consult our wiki:
%   http://control.ee.ethz.ch/~mpt/3/Main/MigrationFromMPT2
%   http://control.ee.ethz.ch/~mpt/3/UI/Control
%   http://control.ee.ethz.ch/~mpt/3/UI/Systems
%   http://control.ee.ethz.ch/~mpt/3/UI/ClosedLoop
%
% Also consult "help MPCController" and "help EMPCController".

mpt_obsoleteFunction;
fprintf('\n');
help(mfilename);

% import model from sysStruct/probStruct
model = mpt_import(sysStruct, probStruct);

% all MPC controllers in MPT3 are on-line by default
ctrl = MPCController(model, probStruct.N);

if nargin==2 || ( nargin==3 && isequal(type, 'explicit') )
	% return an explicit representation
	ctrl = ctrl.toExplicit();
	
	% remove overlaps if possible
	if numel(ctrl.optimizer)>1 && probStruct.norm~=2
		fprintf('Removing overlaps:\n');
		ctrl.optimizer = ctrl.optimizer.min('obj');
	end
end

end
