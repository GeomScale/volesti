#!/usr/bin/octave -qf
% sampling_demo.m
%
% Demonstrates sampling from polytopes with visualization

addpath('../inst');

fprintf('\n=== Sampling Demonstration ===\n\n');

%% Example 1: Sampling from 2D square
fprintf('Example 1: Sampling from 2D Unit Square\n');
fprintf('------------------------------------------\n\n');

% Define unit square [-1, 1]^2
A = [1 0; -1 0; 0 1; 0 -1];
b = ones(4, 1);

% Sample 500 points with different parameters
samples1 = sample_points(A, b, 500, 0, 1, 0);    % CDHR, no burn-in
samples2 = sample_points(A, b, 500, 0, 10, 100);  % CDHR, walk_length=10, burn-in=100

fprintf('Generated %d samples from 2D square\n', size(samples1, 2));
fprintf('Sample statistics:\n');
fprintf('  Mean: [%.3f, %.3f] (expected [0, 0])\n', mean(samples1, 2));
fprintf('  Std:  [%.3f, %.3f] (expected [%.3f, %.3f])\n', ...
        std(samples1, 0, 2), 2/sqrt(12), 2/sqrt(12));

fprintf('\nTo visualize (requires graphics):\n');
fprintf('  figure; plot(samples1(1,:), samples1(2,:), ''.''); axis equal;\n');
fprintf('  title(''Uniform Samples from Unit Square'');\n\n');

%% Example 2: Sampling from 3D simplex
fprintf('Example 2: Sampling from 3D Simplex\n');
fprintf('------------------------------------\n\n');

% Define 3D simplex as V-polytope
V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

% Sample using RDHR (better for V-polytopes)
samples = sample_points(V, 300, 1, 5, 50);

fprintf('Generated %d samples from 3D simplex\n', size(samples, 2));
fprintf('Sample statistics:\n');
fprintf('  Mean: [%.3f, %.3f, %.3f] (expected [0.25, 0.25, 0.25])\n', ...
        mean(samples, 2));
fprintf('  All non-negative: %s\n', all(all(samples >= -1e-10)) ? 'Yes' : 'No');
fprintf('  Sum <= 1: %s\n', all(sum(samples) <= 1 + 1e-10) ? 'Yes' : 'No');

fprintf('\nTo visualize (requires graphics):\n');
fprintf('  figure; plot3(samples(1,:), samples(2,:), samples(3,:), ''.'');\n');
fprintf('  title(''Uniform Samples from 3D Simplex'');\n\n');

%% Example 3: Comparing walk types
fprintf('Example 3: Comparing Different Walk Types\n');
fprintf('------------------------------------------\n\n');

% 3D cube
A = [eye(3); -eye(3)];
b = ones(6, 1);

fprintf('Sampling 200 points from 3D cube with each walk type:\n\n');

% CDHR
tic;
samples_cdhr = sample_points(A, b, 200, 0, 1, 0, false);
time_cdhr = toc;
fprintf('  CDHR:      %.3f seconds\n', time_cdhr);

% RDHR  
tic;
samples_rdhr = sample_points(A, b, 200, 1, 1, 0, false);
time_rdhr = toc;
fprintf('  RDHR:      %.3f seconds\n', time_rdhr);

% Ball Walk
tic;
samples_ball = sample_points(A, b, 200, 2, 1, 0, false);
time_ball = toc;
fprintf('  Ball Walk: %.3f seconds\n\n', time_ball);

fprintf('All walk types should produce similar distributions\n');
fprintf('for symmetric polytopes like the cube.\n\n');

%% Example 4: Effect of walk_length
fprintf('Example 4: Effect of walk_length Parameter\n');
fprintf('-------------------------------------------\n\n');

fprintf('Sampling with different walk_length values:\n\n');

A = [1 0; -1 0; 0 1; 0 -1];
b = ones(4, 1);

for wl = [1, 5, 10, 20]
    tic;
    samples = sample_points(A, b, 100, 0, wl, 0, false);
    t = toc;
    fprintf('  walk_length=%2d: %.3f sec, mean=[%.3f, %.3f]\n', ...
            wl, t, mean(samples, 2));
end

fprintf('\nLarger walk_length:\n');
fprintf('  - Better mixing (more independent samples)\n');
fprintf('  - Higher computation time\n');
fprintf('  - Recommended: 5-10 for most applications\n\n');

%% Example 5: Burn-in period
fprintf('Example 5: Importance of Burn-in Period\n');
fprintf('----------------------------------------\n\n');

fprintf('Starting from corner (non-typical point):\n\n');

% The first few samples might be biased without burn-in
samples_no_burn = sample_points(A, b, 100, 0, 1, 0, false);
samples_burn50 = sample_points(A, b, 100, 0, 1, 50, false);
samples_burn200 = sample_points(A, b, 100, 0, 1, 200, false);

fprintf('  No burn-in:  mean=[%.3f, %.3f]\n', mean(samples_no_burn, 2));
fprintf('  Burn-in=50:  mean=[%.3f, %.3f]\n', mean(samples_burn50, 2));
fprintf('  Burn-in=200: mean=[%.3f, %.3f]\n\n', mean(samples_burn200, 2));

fprintf('Burn-in helps reach equilibrium distribution.\n');
fprintf('Recommended: 10-100 iterations for simple polytopes\n\n');

fprintf('=== Demonstration Complete ===\n\n');
fprintf('For more examples, see the test/ directory.\n');
fprintf('For API documentation, see README.md\n\n');
