#!/usr/bin/octave -qf
% test_volesti.m
%
% Verification script for Volesti Octave interface
% Tests volume computation on a 2D hypercube with known analytical result

% Load the compute_volume function from the .oct plugin
autoload('compute_volume', fullfile(pwd, 'volesti_volume.oct'));

fprintf('\\n=== Volesti Octave Interface - Volume Test ===\\n\\n');

% Define a 2D hypercube (square): -1 <= x <= 1, -1 <= y <= 1
% Constraints: x <= 1, -x <= 1, y <= 1, -y <= 1
% Which is: Ax <= b where
%   A = [ 1  0;     (x <= 1)
%        -1  0;     (-x <= 1, or x >= -1)
%         0  1;     (y <= 1)
%         0 -1]     (-y <= 1, or y >= -1)
%   b = [1; 1; 1; 1]

A = [ 1  0;
     -1  0;
      0  1;
      0 -1];

b = ones(4, 1);

fprintf('Testing 2D Hypercube:\n');
fprintf('  Constraints: -1 <= x <= 1, -1 <= y <= 1\n');
fprintf('  Shape: Square with side length 2\n');
fprintf('  Expected Volume: 2 × 2 = 4.0\n\n');

% Compute volume using Volesti wrapper
fprintf('Computing volume using Volesti...\n');
computed_volume = compute_volume(A, b, 0.1, 10);

% Expected analytical result
expected_volume = 4.0;

% Display results
fprintf('\n--- Results ---\n');
fprintf('Computed Volume: %.4f\n', computed_volume);
fprintf('Expected Volume: %.4f\n', expected_volume);

% Check if results match (with tolerance for numerical approximation)
tolerance = 0.5;  % Volume approximation can have some error
error = abs(computed_volume - expected_volume);

fprintf('\nAbsolute Error: %.4f\n', error);

if error < tolerance
    fprintf('\n SUCCESS: Results match!\n');
    fprintf('The Octave wrapper is working correctly.\n\n');
    exit(0);
else
    fprintf('\n FAILURE: Results do not match.\n');
    fprintf('Error exceeds tolerance of %.2f\n\n', tolerance);
    exit(1);
end
