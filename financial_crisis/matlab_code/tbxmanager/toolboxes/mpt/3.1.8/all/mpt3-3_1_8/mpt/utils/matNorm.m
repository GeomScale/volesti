function n = matNorm(A)
% MATNORM Return norm of the rows of the matrix A
%
% n = matNorm(A)
%

n = sqrt(sum(A.*A,2));

end
