function [p,n,z] = inertia(M)
% Returns the matrix inertia of square matrix M.
% USAGE:
%
%      [p,n,z] = inertia(M)
%
% where p: # of positive eigenvalues
%       n: # of negative eigenvalues
%       z: # of zero eigenvalues

assert(size(M,1)==size(M,2),'Matrix must be square');
e = eig(M);

p = sum(e>0);
n = sum(e<0);
z = sum(e==0);
