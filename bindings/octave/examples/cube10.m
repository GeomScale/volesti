#!/usr/bin/octave -qf
% example_10d.m
%
% Example: Computing the volume of a 10-dimensional hypercube
% This demonstrates the Octave interface can handle higher-dimensional problems

% Load the compute_volume function from the .oct plugin
autoload('compute_volume', fullfile(pwd, 'volesti_volume.oct'));

fprintf('\n=== 10D Hypercube Volume Computation ===\n\n');

% Define a 10-dimensional hypercube
% Side length 2 implies volume = 2^10 = 1024
dim = 10;

% H-representation: Ax <= b
% A hypercube has 2*dim constraints (lower and upper bounds for each dimension)
% Constraints: -1 <= x_i <= 1 for i = 1...dim
A = [eye(dim); -eye(dim)];
b = ones(2*dim, 1);

fprintf('Testing %dD Hypercube:\n', dim);
fprintf('  Constraints: -1 <= x_i <= 1 for i = 1...%d\n', dim);
fprintf('  Shape: Hypercube with side length 2 in each dimension\n');
fprintf('  Expected Volume: 2^%d = %.0f\n\n', dim, 2^dim);

% Compute volume using Volesti wrapper
fprintf('Computing volume using Volesti...\n');
fprintf('(This may take a few moments for higher dimensions)\n\n');

tic;
computed_volume = compute_volume(A, b);
elapsed = toc;

% Expected analytical result
expected_volume = 2^dim;

% Display results
fprintf('--- Results ---\n');
fprintf('Computed Volume: %.4f\n', computed_volume);
fprintf('Expected Volume: %.0f\n', expected_volume);
fprintf('Computation Time: %.2f seconds\n', elapsed);

% Calculate relative error
relative_error = abs(computed_volume - expected_volume) / expected_volume * 100;
fprintf('\nRelative Error: %.2f%%\n', relative_error);

% Tolerance for higher-dimensional approximation
tolerance_percent = 20;  % 20% tolerance for stochastic estimation

if relative_error < tolerance_percent
    fprintf('\n SUCCESS: Result is within acceptable tolerance!\n');
    fprintf('The wrapper successfully handles high-dimensional polytopes.\n\n');
else
    fprintf('\n  WARNING: Error exceeds %.0f%% tolerance.\n', tolerance_percent);
    fprintf('Note: Volume approximation becomes more challenging in high dimensions.\n\n');
end
