function [nA,nb] = normalize(A, b)
% NORMALIZE Return normalized rows of the matrix A
%
% n = normalize(A, b)
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

n = sqrt(sum(A.*A,2));
if n>=MPTOPTIONS.abs_tol
    nA = A ./ repmat(n,1,size(A,2));
    if nargin > 1
        nb = b ./ repmat(n,1,size(b,2));
    end
else
    nA = A;
    if nargin > 1
        nb = b;
    end
end

end
