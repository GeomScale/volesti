function P = Hpolytope(A, b)
% HPOLYTOPE Create an H-polytope (H-representation: Ax <= b)
%
%   P = Hpolytope(A, b) creates an H-polytope struct representing the
%   convex polytope defined by the linear inequalities Ax <= b.
%
% Arguments:
%   A - m x n constraint matrix (m constraints, n dimensions)
%   b - m x 1 constraint vector
%
% Returns:
%   P - Struct with fields:
%       .type - 'Hpolytope'
%       .A    - Constraint matrix
%       .b    - Constraint vector
%       .d    - Dimension (n)
%       .m    - Number of constraints (m)
%
% Example:
%   % Create a 2D unit square: -1 <= x,y <= 1
%   A = [1 0; -1 0; 0 1; 0 -1];
%   b = ones(4, 1);
%   P = Hpolytope(A, b);
%   vol = volume(P);
%
% See also: volume, Vpolytope

    % Input validation
    if nargin < 2
        error('Hpolytope: requires 2 arguments (A, b)');
    end
    
    % Check dimensions
    [m, n] = size(A);
    if size(b, 1) ~= m || size(b, 2) ~= 1
        error('Hpolytope: b must be a column vector with %d rows', m);
    end
    
    % Create polytope struct
    P = struct();
    P.type = 'Hpolytope';
    P.A = A;
    P.b = b;
    P.d = n;  % dimension
    P.m = m;  % number of constraints
end
