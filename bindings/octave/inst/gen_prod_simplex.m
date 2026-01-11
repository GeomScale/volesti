function P = gen_prod_simplex(dimension)
% GEN_PROD_SIMPLEX Generate product of two d-dimensional unit simplices
%
%   P = gen_prod_simplex(d) creates a 2d-dimensional polytope that is the
%   product of two d-dimensional unit simplices in H-representation
%
% Arguments:
%   dimension - The dimension of each simplex
%
% Returns:
%   P - H-polytope struct representing the product
%
% Notes:
%   The product of two d-dimensional simplices is a 2d-dimensional polytope
%   Exact volume is (1/d!)^2
%
% Examples:
%   % Product of two 5-dimensional simplices (10D polytope)
%   P = gen_prod_simplex(5);
%   vol = volesti_volume(P);  % Should be ~(1/120)^2 = 1/14400
%
% See also: gen_cube, gen_cross, gen_simplex, Hpolytope

    % Validate input
    if ~isscalar(dimension) || dimension < 1 || dimension ~= floor(dimension)
        error('gen_prod_simplex: dimension must be a positive integer');
    end
    
    % Product of simplices: S_d × S_d
    % Each d-simplex has d+1 constraints:
    %   x_i >= 0 for i=1..d
    %   sum(x_i) <= 1
    %
    % Product has 2d dimensions with 2(d+1) constraints
    
    d = dimension;
    total_dim = 2 * d;
    num_constraints = 2 * (d + 1);
    
    A = zeros(num_constraints, total_dim);
    b = zeros(num_constraints, 1);
    
    % First simplex constraints (first d dimensions)
    % -x_i <= 0  (i.e., x_i >= 0) for i=1..d
    for i = 1:d
        A(i, i) = -1;
        b(i) = 0;
    end
    
    % Sum constraint: x_1 + ... + x_d <= 1
    A(d+1, 1:d) = 1;
    b(d+1) = 1;
    
    % Second simplex constraints (last d dimensions)
    % -x_i <= 0  (i.e., x_i >= 0) for i=d+1..2d
    for i = 1:d
        A(d+1+i, d+i) = -1;
        b(d+1+i) = 0;
    end
    
    % Sum constraint: x_{d+1} + ... + x_{2d} <= 1
    A(2*d+2, (d+1):(2*d)) = 1;
    b(2*d+2) = 1;
    
    P = Hpolytope(A, b);
end
