function P = gen_cross(dimension, representation)
% GEN_CROSS Generate a d-dimensional cross polytope
%
%   P = gen_cross(d, type) creates a d-dimensional cross polytope in H- or V-representation
%
% Arguments:
%   dimension - The dimension of the cross polytope
%   representation - 'H' for H-representation or 'V' for V-representation (default: 'H')
%
% Returns:
%   P - Polytope struct (Hpolytope or Vpolytope)
%
% Notes:
%   Cross polytope is defined by: sum(|x_i|) <= 1
%   In d dimensions, has 2^d facets and 2d vertices
%
% Examples:
%   % 5D cross polytope in H-representation
%   P = gen_cross(5, 'H');
%   
%   % 15D cross polytope in V-representation  
%   P = gen_cross(15, 'V');
%
% See also: gen_cube, gen_simplex, Hpolytope, Vpolytope

    % Default to H-representation
    if nargin < 2
        representation = 'H';
    end
    
    % Validate inputs
    if ~isscalar(dimension) || dimension < 1 || dimension ~= floor(dimension)
        error('gen_cross: dimension must be a positive integer');
    end
    
    if strcmpi(representation, 'H')
        % H-representation: All combinations of ±1 for coefficients
        % sum(±x_i) <= 1 for all 2^d combinations of signs
        
        num_constraints = 2^dimension;
        A = zeros(num_constraints, dimension);
        b = ones(num_constraints, 1);
        
        % Generate all sign combinations
        for i = 0:(num_constraints-1)
            for j = 1:dimension
                % Extract j-th bit and map 0->-1, 1->+1
                bit = bitget(i, j);
                A(i+1, j) = 2*bit - 1;
            end
        end
        
        P = Hpolytope(A, b);
        
    elseif strcmpi(representation, 'V')
        % V-representation: 2d vertices at ±e_i (standard basis vectors)
        
        V = zeros(2 * dimension, dimension);
        
        % Positive and negative unit vectors
        for i = 1:dimension
            V(i, i) = 1;               % +e_i
            V(dimension + i, i) = -1;  % -e_i
        end
        
        P = Vpolytope(V);
        
    else
        error('gen_cross: representation must be ''H'' or ''V''');
    end
end
