#!/usr/bin/octave -qf
% test_generators.m
%
% Test polytope generators with known volumes

fprintf('\n=== Testing Polytope Generators ===\n\n');

% Add inst/ to path
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'inst'));

% Load MEX function
autoload('compute_volume', fullfile(fileparts(mfilename('fullpath')), '..', 'inst', 'volume.oct'));

%% Test 1: 3D Cube (H-representation)
fprintf('Test 1: 3D Cube (H-representation)\n');
P = gen_cube(3, 'H');
fprintf('  Generated: %d constraints, %d dimensions\n', P.m, P.d);
vol = volesti_volume(P, 0.1, 10, false);
expected = 8.0;  % 2^3
error = abs(vol - expected) / expected * 100;
fprintf('  Volume: %.4f (expected: %.4f, error: %.2f%%)\n', vol, expected, error);
if error < 10
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
end

%% Test 2: 3D Cross (H-representation)
fprintf('Test 2: 3D Cross Polytope (H-representation)\n');
P = gen_cross(3, 'H');
fprintf('  Generated: %d constraints, %d dimensions\n', P.m, P.d);
vol = volesti_volume(P, 0.1, 10, false);
expected = 2^3 / factorial(3);  % 8/6 = 1.333...
error = abs(vol - expected) / expected * 100;
fprintf('  Volume: %.4f (expected: %.4f, error: %.2f%%)\n', vol, expected, error);
if error < 10
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
end

%% Test 3: Product Simplex (3D)
fprintf('Test 3: Product of 3D Simplices (6D polytope)\n');
P = gen_prod_simplex(3);
fprintf('  Generated: %d constraints, %d dimensions\n', P.m, P.d);
vol = volesti_volume(P, 0.1, 10, false);
expected = (1/factorial(3))^2;  % (1/6)^2 = 0.027778
error = abs(vol - expected) / expected * 100;
fprintf('  Volume: %.6f (expected: %.6f, error: %.2f%%)\n', vol, expected, error);
if error < 15
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
end

%% Test 4: 3D Cube (V-representation)
fprintf('Test 4: 3D Cube (V-representation)\n');
P = gen_cube(3, 'V');
fprintf('  Generated: %d vertices, %d dimensions\n', P.m, P.d);
vol = volesti_volume(P, 0.1, 10, false);
expected = 8.0;
error = abs(vol - expected) / expected * 100;
fprintf('  Volume: %.4f (expected: %.4f, error: %.2f%%)\n', vol, expected, error);
if error < 10
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
end

%% Test 5: 3D Cross (V-representation)
fprintf('Test 5: 3D Cross Polytope (V-representation)\n');
P = gen_cross(3, 'V');
fprintf('  Generated: %d vertices, %d dimensions\n', P.m, P.d);
vol = volesti_volume(P, 0.1, 10, false);
expected = 2^3 / factorial(3);
error = abs(vol - expected) / expected * 100;
fprintf('  Volume: %.4f (expected: %.4f, error: %.2f%%)\n', vol, expected, error);
if error < 10
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
end

fprintf('=== Generator Tests Complete ===\n\n');
