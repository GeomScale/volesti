function obj = Penalty(Q, type)
% Deprecated! Use @Function classes to define penalties.
%
% x'*Q*x:      P = QuadFunction(Q)
% ||Q*x||_1:   P = OneNormFunction(Q)
% ||Q*x||_inf: P = InfNormFunction(Q)
% zero-norm:   P = AffFunction(Q)

if nargin < 2
	type = -1;
end

switch type,
	case 1,
		s = 'Use obj.penalty = OneNormFunction(Q) to define ||Q*z||_1 penalty.';
	case 2,
		s = 'Use obj.penalty = QuadFunction(Q) to define z''*Q*z penalty.';
	case Inf,
		s = 'Use obj.penalty = InfNormFunction(Q) to define ||Q*z||_inf penalty.';
	case 0,
		s = 'Use obj.penalty = AffFunction(Q) to define Q*z penalty.';
	otherwise
		s = 'Use QuadFunction, OneNormFunction or InfNormFunction instead.';
end

fprintf('%s\n', repmat('=', 1, max(20, length(s))));
fprintf('The Penalty object is deprecated.\n');
fprintf('%s\n', s);
fprintf('%s\n', repmat('=', 1, max(20, length(s))));

switch type,
	case 1,
		obj = OneNormFunction(Q);
	case 2,
		obj = QuadFunction(Q);
	case Inf,
		obj = InfNormFunction(Q);
	case 0,
		% just for backwards compatibility
		obj = AffFunction(Q);
	case -1
		error('Cannot continue.');
	otherwise,
		error('Unrecognized norm type "%s".', num2str(type));
end

end
