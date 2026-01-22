function P = gen_simplex(dimension, representation)
%GEN_SIMPLEX Generate d-dimensional unit simplex
%
% P = gen_simplex(dimension, representation)
%
% Generates the d-dimensional unit simplex in H- or V-representation.
% The unit simplex is defined as {x : x_i >= 0, sum(x_i) <= 1}
%
% Parameters:
%   dimension : The dimension of the simplex
%   representation : (optional) 'H' for H-representation (default)
%                               'V' for V-representation
%
% Returns:
%   P : Polytope struct (Hpolytope or Vpolytope)
%       Volume: 1/d! (volume of unit simplex)
%
% Examples:
%   % 10-dimensional simplex in H-representation
%   P = gen_simplex(10, 'H');
%
%   % 20-dimensional simplex in V-representation
%   P = gen_simplex(20, 'V');
%
%   % Default is H-representation
%   P = gen_simplex(5);

    % Default representation is H
    if nargin < 2
        representation = 'H';
    end
    
    % Validate inputs
    if ~isnumeric(dimension) || dimension < 1 || floor(dimension) ~= dimension
        error('gen_simplex: dimension must be a positive integer');
    end
    
    if ~ischar(representation) || ~(strcmp(representation, 'H') || strcmp(representation, 'V'))
        error('gen_simplex: representation must be ''H'' or ''V''');
    end
    
    % Compute volume (1 / d!)
    volume = 1 / factorial(dimension);
    
    if strcmp(representation, 'V')
        % V-representation: vertices are origin + standard basis vectors
        % Vertices: (0,0,...,0), (1,0,...,0), (0,1,0,...,0), ..., (0,0,...,1)
        V = [zeros(1, dimension); eye(dimension)];
        P = Vpolytope(V);
        P.volume = volume;
    else
        % H-representation: {x : x_i >= 0, sum(x_i) <= 1}
        % Constraints:
        %   -x_1 <= 0  (x_1 >= 0)
        %   -x_2 <= 0  (x_2 >= 0)
        %   ...
        %   x_1 + x_2 + ... + x_d <= 1
        
        A = [-eye(dimension); ones(1, dimension)];
        b = [zeros(dimension, 1); 1];
        P = Hpolytope(A, b);
        P.volume = volume;
    end
end
