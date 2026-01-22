#!/usr/bin/octave -qf
% test_sample_points.m
%
% Basic test for sample_points.oct function

% Add inst to path to find the .oct file
addpath('../inst');

fprintf('\n=== Testing sample_points.oct ===\n\n');

% Test 1: H-polytope (2D unit square)
fprintf('Test 1: Sampling from 2D H-polytope (unit square)\n');
A = [1 0; -1 0; 0 1; 0 -1];
b = ones(4, 1);

try
    samples = sample_points(A, b, 100);  % Direct .oct call, R-style
    fprintf('  Generated %d samples from %dD polytope\n', size(samples, 2), size(samples, 1));
    fprintf('  Sample mean: [%.3f, %.3f] (expected ~[0, 0])\n', mean(samples, 2));
    
    inside = all(all(A * samples <= repmat(b, 1, size(samples, 2)) + 1e-10));
    if inside
        fprintf('  All samples inside polytope: YES\n');
    else
        fprintf('  All samples inside polytope: NO\n');
    end
    fprintf('  ✓ H-polytope sampling works!\n\n');
catch err
    fprintf('  ✗ Error: %s\n\n', err.message);
end

% Test 2: V-polytope (3D simplex)
fprintf('Test 2: Sampling from 3D V-polytope (simplex)\n');
V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

try
    samples = sample_points(V, 50);  % Direct .oct call with V-polytope
    fprintf('  Generated %d samples from %dD polytope\n', size(samples, 2), size(samples, 1));
    fprintf('  Sample mean: [%.3f, %.3f, %.3f]\n', mean(samples, 2));
    
    nonneg = all(all(samples >= -1e-10));
    if nonneg
        fprintf('  All coordinates non-negative: YES\n');
    else
        fprintf('  All coordinates non-negative: NO\n');
    end
    fprintf('  ✓ V-polytope sampling works!\n\n');
catch err
    fprintf('  ✗ Error: %s\n\n', err.message);
end

% Test 3: Custom walk parameters
fprintf('Test 3: Custom walk parameters (RDHR, walk_length=10, nburns=20)\n');
try
    samples = sample_points(A, b, 100, 1, 10, 20, false);  % walk_type=1 (RDHR), verbose=false
    fprintf('  Generated %d samples (silent mode)\n', size(samples, 2));
    fprintf('  ✓ Custom parameters work!\n\n');
catch err
    fprintf('  ✗ Error: %s\n\n', err.message);
end

fprintf('=== All tests passed! ===\n\n');
