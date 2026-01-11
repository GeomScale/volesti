#!/usr/bin/octave -qf
% test_Hvol.m
%
% Volume test for H-polytopes (matches Rvolesti's test_Hvol.R)
% Tests volume computation on a 2D hypercube with known analytical result

fprintf('\n=== Octave-Volesti: H-Polytope Volume Test ===\n\n');

% Add inst/ to path for API functions
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'inst'));

% Load the compiled volume function
autoload('compute_volume', fullfile(fileparts(mfilename('fullpath')), '..', 'inst', 'volume.oct'));

% Define a 2D hypercube (square): -1 <= x <= 1, -1 <= y <= 1
% Constraints: x <= 1, -x <= 1, y <= 1, -y <= 1
A = [ 1  0;
     -1  0;
      0  1;
      0 -1];

b = ones(4, 1);

fprintf('Testing 2D Hypercube:\n');
fprintf('  Constraints: -1 <= x <= 1, -1 <= y <= 1\n');
fprintf('  Shape: Square with side length 2\n');
fprintf('  Expected Volume: 2 × 2 = 4.0\n\n');

% Create H-polytope using new API (matches Rvolesti)
P = Hpolytope(A, b);
fprintf('Created H-polytope: %d constraints, %d dimensions\n', P.m, P.d);

% Compute volume using new API (matches R: volume(P))
fprintf('\nComputing volume using Volesti SOB algorithm...\n');
computed_volume = volesti_volume(P, 0.1, 10);

% Expected analytical result
expected_volume = 4.0;

% Display results
fprintf('\n--- Results ---\n');
fprintf('Computed Volume: %.4f\n', computed_volume);
fprintf('Expected Volume: %.4f\n', expected_volume);

% Check if results match (with tolerance for numerical approximation)
tolerance = 0.5;  % Volume approximation can have some error
error = abs(computed_volume - expected_volume);
rel_error = error / expected_volume;

fprintf('\nAbsolute Error: %.4f\n', error);
fprintf('Relative Error: %.2f%%\n', rel_error * 100);

if error < tolerance
    fprintf('\n✓ SUCCESS: Results match!\n');
    fprintf('The Octave wrapper is working correctly.\n\n');
    exit(0);
else
    fprintf('\n✗ FAILURE: Results do not match.\n');
    fprintf('Error exceeds tolerance of %.2f\n\n', tolerance);
    exit(1);
end
