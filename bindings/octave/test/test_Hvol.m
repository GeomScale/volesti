#!/usr/bin/octave -qf
% test_Hvol.m
%
% H-polytope volume test matching Rvolesti's test_Hvol.R
% Tests volume computation on standard polytopes with known exact volumes

fprintf('\n=== H-Polytope Volume Test (matching Rvolesti) ===\n\n');

% Add inst/ to path
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'inst'));

% Load MEX function
autoload('compute_volume', fullfile(fileparts(mfilename('fullpath')), '..', 'inst', 'volume.oct'));

% Test parameters matching Rvolesti
num_of_exps = 5;  % Reduced from 10 for faster testing
tolerance = 0.2;   % 20% tolerance (same as Rvolesti)
seed = 5;          % Fixed seed for reproducibility
algo = 'SOB';      % We only have SOB algorithm

fprintf('Algorithm: %s\n', algo);
fprintf('Number of experiments: %d\n', num_of_exps);
fprintf('Tolerance: %.1f%%\n\n', tolerance * 100);

all_passed = true;

%% Test 1: 10D Hypercube
fprintf('Test 1: Volume H-cube(10)\n');
P = gen_cube(10, 'H');
fprintf('  Polytope: %d constraints, %d dimensions\n', P.m, P.d);

exactvol = 2^10;  % 1024
vol = 0;
for j = 1:num_of_exps
    vol = vol + volesti_volume(P, 1.0, 1, false);  % SOB defaults
end
vol = vol / num_of_exps;
error = abs(vol - exactvol) / exactvol;

fprintf('  Computed: %.2f, Expected: %.2f\n', vol, exactvol);
fprintf('  Relative error: %.2f%%\n', error * 100);

if error < tolerance
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
    all_passed = false;
end

%% Test 2: 5D Cross Polytope
fprintf('Test 2: Volume H-cross(5)\n');
P = gen_cross(5, 'H');
fprintf('  Polytope: %d constraints, %d dimensions\n', P.m, P.d);

exactvol = 2^5 / factorial(5);  % 32/120 = 0.2666667
vol = 0;
for j = 1:num_of_exps
    vol = vol + volesti_volume(P, 1.0, 1, false);
end
vol = vol / num_of_exps;
error = abs(vol - exactvol) / exactvol;

fprintf('  Computed: %.4f, Expected: %.4f\n', vol, exactvol);
fprintf('  Relative error: %.2f%%\n', error * 100);

if error < tolerance
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
    all_passed = false;
end

%% Test 3: Product of 5D Simplices
fprintf('Test 3: Volume H-prod_simplex(5,5)\n');
P = gen_prod_simplex(5);
fprintf('  Polytope: %d constraints, %d dimensions\n', P.m, P.d);

exactvol = (1/factorial(5))^2;  % (1/120)^2
vol = 0;
for j = 1:num_of_exps
    vol = vol + volesti_volume(P, 1.0, 1, false);
end
vol = vol / num_of_exps;
error = abs(vol - exactvol) / exactvol;

fprintf('  Computed: %.8f, Expected: %.8f\n', vol, exactvol);
fprintf('  Relative error: %.2f%%\n', error * 100);

if error < tolerance
    fprintf('  ✓ PASS\n\n');
else
    fprintf('  ✗ FAIL\n\n');
    all_passed = false;
end

%% Summary
fprintf('=== Test Summary ===\n');
if all_passed
    fprintf('All tests PASSED! ✓\n\n');
    exit(0);
else
    fprintf('Some tests FAILED! ✗\n\n');
    exit(1);
end
