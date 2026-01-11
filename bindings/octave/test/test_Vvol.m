#!/usr/bin/octave -qf
% test_Vvol.m
%
% Volume test for V-polytopes (matches Rvolesti's test_Vvol.R)
% Tests volume computation on a 3D cube defined by 8 vertices

fprintf('\n=== Octave-Volesti: V-Polytope Volume Test ===\n\n');

% Add inst/ to path for API functions
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'inst'));

% Load the compiled volume function
volesti_oct_path = fullfile(fileparts(mfilename('fullpath')), '..', 'inst', 'volume.oct');
autoload('compute_volume', volesti_oct_path);

% Define a 3D unit cube by its 8 vertices
% Cube: [-1, 1]^3
fprintf('Testing 3D Cube (V-representation):\n');
fprintf('  8 vertices defining cube [-1,1]^3\n');
fprintf('  Expected Volume: 2 × 2 × 2 = 8.0\n\n');

% Create vertex matrix (each row is a vertex)
V = [ 1  1  1;   %  (+,+,+)
      1  1 -1;   %  (+,+,-)
      1 -1  1;   %  (+,-,+)
      1 -1 -1;   %  (+,-,-)
     -1  1  1;   %  (-,+,+)
     -1  1 -1;   %  (-,+,-)
     -1 -1  1;   %  (-,-,+)
     -1 -1 -1 ]; %  (-,-,-)

% Create V-polytope using new API (matches Rvolesti)
P = Vpolytope(V);
fprintf('Created V-polytope: %d vertices, %d dimensions\n', P.m, P.d);

% Compute volume using new API
fprintf('\nComputing volume using Volesti SOB algorithm...\n');
computed_volume = volesti_volume(P, 0.1, 10);

% Expected analytical result
expected_volume = 8.0;

% Display results
fprintf('\n--- Results ---\n');
fprintf('Computed Volume: %.4f\n', computed_volume);
fprintf('Expected Volume: %.4f\n', expected_volume);

% Check if results match (with tolerance for numerical approximation)
tolerance = 1.0;  % V-polytope volume can have larger error
error = abs(computed_volume - expected_volume);
rel_error = error / expected_volume;

fprintf('\nAbsolute Error: %.4f\n', error);
fprintf('Relative Error: %.2f%%\n', rel_error * 100);

if error < tolerance
    fprintf('\n✓ SUCCESS: Results match!\n');
    fprintf('The V-polytope wrapper is working correctly.\n\n');
    exit(0);
else
    fprintf('\n✗ FAILURE: Results do not match.\n');
    fprintf('Error exceeds tolerance of %.2f\n\n', tolerance);
    exit(1);
end
