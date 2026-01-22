function P = gen_birkhoff(n)
%GEN_BIRKHOFF Generate n-Birkhoff polytope
%
% P = gen_birkhoff(n)
%
% Generates the n-Birkhoff polytope in H-representation.
% The Birkhoff polytope B_n is the set of n×n doubly stochastic matrices.
% Dimension: (n-1)^2 (full dimensional in flattened space)
%
% A doubly stochastic matrix satisfies:
%   - All elements >= 0
%   - Each row sums to 1
%   - Each column sums to 1
%
% Parameters:
%   n : The order of the Birkhoff polytope (must be >= 2)
%
% Returns:
%   P : H-polytope struct representing the Birkhoff polytope
%      

% Dimension: (n-1)^2
%
% Examples:
%   % Generate the 5-Birkhoff polytope (dimension 16)
%   P = gen_birkhoff(5);
%
%   % Generate the 3-Birkhoff polytope (dimension 4)
%   P = gen_birkhoff(3);
%
% Note:
%   - Vertices are permutation matrices
%   - Volume is known for small n
%   - Only H-representation is supported

    % Validate input
    if ~isnumeric(n) || n < 2 || floor(n) ~= n
        error('gen_birkhoff: n must be an integer >= 2');
    end
    
    % For now, use a simple implementation for n=2,3
    % For larger n, would need C++ generator
    
    if n == 2
        % B_2 is a single point: [[0.5, 0.5], [0.5, 0.5]]
        % In 1D: just the value 0.5
        A = [1; -1];
        b = [0.5; -0.5];
    elseif n == 3
        % B_3 is 4-dimensional
        % Constraints are more complex, use proper formulation
        % For a 3×3 doubly stochastic matrix:
        % [x11 x12 x13]
        % [x21 x22 x23]
        % [x31 x32 x33]
        %
        % Parameterized by (x11, x12, x21, x22):
        % x13 = 1 - x11 - x12
        % x23 = 1 - x21 - x22  
        % x31 = 1 - x11 - x21
        % x32 = 1 - x12 - x22
        % x33 = 1 - x13 - x23 = 1 - (1-x11-x12) - (1-x21-x22) = x11 + x12 + x21 + x22 - 1
        
        % Constraints:
        % x11, x12, x21, x22 >= 0
        % x11 + x12 <= 1 (so x13 >= 0)
        % x21 + x22 <= 1 (so x23 >= 0)
        % x11 + x21 <= 1 (so x31 >= 0)
        % x12 + x22 <= 1 (so x32 >= 0)
        % x11 + x12 + x21 + x22 >= 1 (so x33 >= 0)
        
        A = [
            -1  0  0  0;  % -x11 <= 0
             0 -1  0  0;  % -x12 <= 0
             0  0 -1  0;  % -x21 <= 0
             0  0  0 -1;  % -x22 <= 0
             1  1  0  0;  % x11 + x12 <= 1
             0  0  1  1;  % x21 + x22 <= 1
             1  0  1  0;  % x11 + x21 <= 1
             0  1  0  1;  % x12 + x22 <= 1
            -1 -1 -1 -1;  % -(x11+x12+x21+x22) <= -1 (x33 >= 0)
        ];
        b = [0; 0; 0; 0; 1; 1; 1; 1; -1];
    else
        % For larger n, this becomes complex
        % Note: In practice, should use C++ generator from volesti
        warning('gen_birkhoff: Full implementation for n>3 requires C++ generator');
        error('gen_birkhoff: Currently only n=2 and n=3 are fully supported');
    end
    
    P = Hpolytope(A, b);
    % Volume formula for Birkhoff polytopes is complex, don't set it
end
