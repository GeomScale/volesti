function P = Vpolytope(V)
% VPOLYTOPE Create a V-polytope (V-representation: convex hull of vertices)
%
%   P = Vpolytope(V) creates a V-polytope struct representing the convex hull
%   of the given vertices.
%
% Arguments:
%   V - num_vertices x n matrix where each row is a vertex in n-dimensional space
%
% Returns:
%   P - Struct with fields:
%       .type - 'Vpolytope'
%       .V    - Vertex matrix (num_vertices x n)
%       .d    - Dimension (n)
%       .m    - Number of vertices
%
% Example:
%   % Create a 2D triangle with vertices at (0,0), (1,0), (0,1)
%   V = [0 0; 1 0; 0 1];
%   P = Vpolytope(V);
%   vol = volesti_volume(P);
%
%   % Create a 3D cube with 8 vertices
%   V = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1;
%        -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
%   P = Vpolytope(V);
%   vol = volesti_volume(P);  % Should be ~8
%
% See also: volesti_volume, Hpolytope

    % Input validation
    if nargin < 1
        error('Vpolytope: requires 1 argument (V)');
    end
    
    % Check that V is a 2D matrix
    if ~ismatrix(V) || ndims(V) ~= 2
        error('Vpolytope: V must be a 2D matrix');
    end
    
    % Get dimensions
    [m, n] = size(V);
    
    if m < n + 1
        warning('Vpolytope: Need at least %d vertices for %dD polytope (got %d)', n+1, n, m);
    end
    
    % Create polytope struct
    P = struct();
    P.type = 'Vpolytope';
    P.V = V;
    P.d = n;  % dimension
    P.m = m;  % number of vertices
end
