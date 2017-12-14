function cost = buildCost(obj, f, hess)
% Build a cost function to match the optMat structure
% Convert min x'*hess*x + f'*x
% to a form to match the optMat structure

if obj.optMat.type == 'V'
    cost.f = [f;zeros(size(obj.optMat.Ae,2)-obj.Dim,1)];
    if nargin > 2
        cost.H = blkdiag(hess, zeros(size(obj.optMat.Ae,2)-obj.Dim));
    end
else
    cost.f = f;
    if nargin > 2,
        cost.H = hess;
    end
end

end
