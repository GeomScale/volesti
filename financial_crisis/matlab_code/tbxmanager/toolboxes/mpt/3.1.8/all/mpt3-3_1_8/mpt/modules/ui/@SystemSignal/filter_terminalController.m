function out = terminalSet(obj, from, args)
% Adds terminal set/penalty

error('Not yet implemented.');

% determine required number of input arguments
if nargin==1
	out = obj.emptyFilterArguments();
	out.parser.addRequired('controller', @(x) isa(x, 'EMPCController'));
	out.incompatible_with = {'terminalPenalty', 'terminalSet'};
	out.execution = 'immediate';
    return
end


% add terminal set
obj.addFilter('terminalSet', args.controller.optimizer.Set);

% add terminal penalty
obj.addFilter('terminalPenalty', args.controller.optimizer.getFunction('obj'));
