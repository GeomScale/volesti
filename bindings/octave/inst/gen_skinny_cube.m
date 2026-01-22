function P = gen_skinny_cube(dimension)
%GEN_SKINNY_CUBE Generate d-dimensional skinny hypercube
%
% P = gen_skinny_cube(dimension)
%
% Generates a d-dimensional skinny hypercube: [-1,1]^(d-1) × [-100,100]
% This is a regular cube in the first (d-1) dimensions and stretched 100x
% in the last dimension, useful for testing rounding algorithms.
%
% Parameters:
%   dimension : The dimension of the skinny hypercube (must be >= 2)
%
% Returns:
%   P : H-polytope struct representing the skinny hypercube
%       Volume: 2^(d-1) * 200
%
% Examples:
%   % Generate a 10-dimensional skinny hypercube
%   P = gen_skinny_cube(10);
%
%   % Generate a 5-dimensional skinny cube
%   P = gen_skinny_cube(5);
%
% Note:
%   Only H-representation is supported (skinny cubes are specified via
%   constraints). The polytope is always full-dimensional.

    % Validate input
    if ~isnumeric(dimension) || dimension < 2 || floor(dimension) ~= dimension
        error('gen_skinny_cube: dimension must be an integer >= 2');
    end
    
    % Create constraints for skinny hypercube
    % First (d-1) dimensions: -1 <= x_i <= 1
    % Last dimension: -100 <= x_d <= 100
    
    % Standard constraints for first (d-1) dimensions
    A_regular = [eye(dimension-1); -eye(dimension-1)];
    b_regular = ones(2*(dimension-1), 1);
    
    % Constraints for last dimension (stretched)
    A_skinny = [zeros(2, dimension-1), [1; -1]];
    b_skinny = [100; 100];
    
    % Combine constraints
    A = [A_regular, zeros(2*(dimension-1), 1); A_skinny];
    b = [b_regular; b_skinny];
    
    % Compute volume: 2^(d-1) * 200
    volume = 2^(dimension-1) * 200;
    
    P = Hpolytope(A, b);
    P.volume = volume;
end
