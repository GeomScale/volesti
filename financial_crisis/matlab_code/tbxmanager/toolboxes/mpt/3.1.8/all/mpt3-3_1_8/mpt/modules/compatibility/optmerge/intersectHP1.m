function m = intersectHP1(H, K, a, b, lpsolver, tol);
% Tobias Geyer, 2003-2004
% Does the hyperplane a*x=b intersect with the polyhedron H*x<=K?
%
% inputs:
% H, K: polyhedron: H*x <= K
% a, b: hyperplane: a*x = b, where the norm of a is equal to 1
% lpsolver
% tolerance
% 
% output:
% m = 0: the hyperplane intersects with the polyhedron
%    +1: the polyhedron lies on the + side of the hyperplane (and does not intersect)
%    -1: the polyhedron lies on the - side of the hyperplane (and does not intersect)
%
% remarks:
% * The tolerance is used to relax the (distances in the) inequalities,
%   i.e. we want to avoid, that the algorithm determines that the hyperplane 
%   cuts the polyhedron, when the hyperplane is only touching it.
% * Relying on the fact, that the hyperplane a*x=b is normed allows us to
%   calculate the distance to the Chebycheff center (this is new compared to ver. 0)
%   speeding up the function

% default tolerance is 0
if nargin < 6
    tol = 0;
end;

[center, R] = polyinnerball(H, K, lpsolver);

if R <= -tol, error('polyhedron is infeasible'); end;

if abs(norm(a,2)-1) > 1e-6, 
    b = b/norm(a, 2); a = a./norm(a, 2);
%     error('hyperplane is not normed'); 
end;

if abs(a*center-b) < R
    % the hyperplane cuts the Chebycheff ball and thus the polyhedron 
    % (what are we lucky... we have just needed one LP to determine this)
    m = 0;  %disp('cuts ball')

elseif a*center > b
    % the center is on the '+' side of the hyperplane
    % so minimize along the normal vector of the hyperplane
    %      min a'x 
    % s.t. Hx <= K
    [xopt, lambda, how] = mpt_solveLP(a, H, K, [], [], [], lpsolver);
    %if ~strcmp(how, 'ok'), warning('LP could not be solved'); end;
    
    if a*xopt >= b-tol  % is this correct?
        % the polyhedron lies completely on the '+' side of the hyperplane
        m = +1; %disp('+')
    else
        % the hyperplane cuts the polyhedron into two halfs
        m = 0;  %disp('cuts')
    end;
      
else
    % the center is on the '-' side of the hyperplane
    % so maximize along the normal vector of the hyperplane
    %      max a'x 
    % s.t. Hx <= K
    [xopt, lambda, how] = mpt_solveLP(-a, H, K, [], [], [], lpsolver);
    %if how ~= 'ok', warning('LP could not be solved'); end;
    %if ~strcmp(how, 'ok'), warning('LP could not be solved'); end;
    
    if a*xopt <= b+tol  % this should be correct...
        % the polyhedron lies completely on the '-' side of the hyperplane
        m = -1; %disp('-')
    else
        % the hyperplane cuts the polyhedron into two halfs
        m = 0;  %disp('cuts')
    end;
end;
