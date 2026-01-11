function result = volesti_volume(P, varargin)
% VOLESTI_VOLUME Compute the volume of a polytope
%
%   vol = volesti_volume(P) computes the volume of polytope P using default parameters
%   vol = volesti_volume(P, epsilon) sets error tolerance (default: 1.0)
%   vol = volesti_volume(P, epsilon, walk_length) sets random walk length (default: 1)
%   vol = volesti_volume(P, epsilon, walk_length, verbose) enables/disables output
%
% Arguments:
%   P - Polytope struct (created with Hpolytope, Vpolytope, etc.)
%   epsilon - (optional) Error tolerance (default: 1.0)
%             Smaller values = more accurate but slower
%   walk_length - (optional) Random walk length (default: 1)
%                 Larger values = better mixing but slower
%   verbose - (optional) Show progress messages (default: true)
%
% Returns:
%   vol - Estimated volume (scalar)
%
% Examples:
%   % Create and compute volume of 3D cube
%   A = [eye(3); -eye(3)];
%   b = ones(6, 1);
%   P = Hpolytope(A, b);
%   vol = volesti_volume(P);  % Should be ~8
%
%   % High accuracy computation
%   vol = volesti_volume(P, 0.1, 10);
%
% See also: Hpolytope, Vpolytope

    % Check that P is a valid polytope struct
    if ~isstruct(P) || ~isfield(P, 'type')
        error('volesti_volume: First argument must be a polytope struct');
    end
    
    % Dispatch based on polytope type
    switch P.type
        case 'Hpolytope'
            result = volume_hpolytope(P, varargin{:});
        case 'Vpolytope'
            result = volume_vpolytope(P, varargin{:});
        case 'Zonotope'
            error('volesti_volume: Zonotope support not yet implemented');
        otherwise
            error('volesti_volume: Unknown polytope type "%s"', P.type);
    end
end

function vol = volume_hpolytope(P, varargin)
    % Internal function to compute H-polytope volume
    
    % Extract polytope data
    A = P.A;
    b = P.b;
    
    % Parse optional arguments
    epsilon = 1.0;
    walk_length = 1;
    verbose = true;
    
    if nargin >= 2
        epsilon = varargin{1};
    end
    if nargin >= 3
        walk_length = varargin{2};
    end
    if nargin >= 4
        verbose = varargin{3};
    end
    
    % Call the compiled MEX function
    % Note: compute_volume is the internal MEX function name
    vol = compute_volume(A, b, epsilon, walk_length, verbose);
end

function vol = volume_vpolytope(P, varargin)
    % Internal function to compute V-polytope volume
    
    % Extract polytope data
    V = P.V;
    
    % Parse optional arguments
    epsilon = 1.0;
    walk_length = 1;
    verbose = true;
    
    if nargin >= 2
        epsilon = varargin{1};
    end
    if nargin >= 3
        walk_length = varargin{2};
    end
    if nargin >= 4
        verbose = varargin{3};
    end
    
    % Call the compiled MEX function with V-polytope signature
    % Pass only V matrix (not A, b)
    vol = compute_volume(V, epsilon, walk_length, verbose);
end
