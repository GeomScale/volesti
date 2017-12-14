function [status, ret] = fast_isEmptySet(H, He)
% returns true if {x | H*[x; -1]<=0, He*[x; -1]=0} is empty

global MPTOPTIONS

S.A = H(:, 1:end-1);
S.b = H(:, end);
if nargin==2
	S.Ae = He(:, 1:end-1);
	S.be = He(:, end);
else
	S.Ae = []; S.be = [];
end
S.lb = []; S.ub = [];
S.f = zeros(size(H, 2)-1, 1); % feasibility problem
S.quicklp = true;
ret = mpt_solve(S);

status = ~(ret.exitflag == MPTOPTIONS.OK);

end
