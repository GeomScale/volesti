function P = gen_cube(dimension, representation)
% GEN_CUBE Generate a d-dimensional unit hypercube [-1,1]^d
%
%   P = gen_cube(d, type) creates a d-dimensional hypercube in H- or V-representation
%
% Arguments:
%   dimension - The dimension of the hypercube
%   representation - 'H' for H-representation or 'V' for V-representation (default: 'H')
%
% Returns:
%   P - Polytope struct (Hpolytope or Vpolytope)
%
% Examples:
%   % 10D hypercube in H-representation
%   P = gen_cube(10, 'H');
%   
%   % 5D hypercube in V-representation
%   P = gen_cube(5, 'V');
%
% See also: gen_cross, gen_simplex, Hpolytope, Vpolytope

    % Default to H-representation
    if nargin < 2
        representation = 'H';
    end
    
    % Validate inputs
    if ~isscalar(dimension) || dimension < 1 || dimension ~= floor(dimension)
        error('gen_cube: dimension must be a positive integer');
    end
    
    if strcmpi(representation, 'H')
        % H-representation: Ax <= b where cube is [-1,1]^d
        % Constraints: x_i <= 1 and -x_i <= 1 for each dimension
        % A = [I; -I], b = ones(2d, 1)
        
        A = [eye(dimension); -eye(dimension)];
        b = ones(2 * dimension, 1);
        
        P = Hpolytope(A, b);
        
    elseif strcmpi(representation, 'V')
        % V-representation: 2^d vertices at all corners {-1,1}^d
        
        num_vertices = 2^dimension;
        V = zeros(num_vertices, dimension);
        
        % Generate all binary combinations
        for i = 0:(num_vertices-1)
            for j = 1:dimension
                % Extract j-th bit and map 0->-1, 1->1
                bit = bitget(i, j);
                V(i+1, j) = 2*bit - 1;
            end
        end
        
        P = Vpolytope(V);
        
    else
        error('gen_cube: representation must be ''H'' or ''V''');
    end
end
