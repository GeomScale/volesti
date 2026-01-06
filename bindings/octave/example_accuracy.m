#!/usr/bin/octave -qf
% example_accuracy.m
%
% Demonstrates the accuracy vs speed tradeoff using optional parameters

% Load the compute_volume function from the .oct plugin
script_dir = fileparts(mfilename('fullpath'));
autoload('compute_volume', fullfile(script_dir, 'volesti_volume.oct'));
fprintf('\n=== Accuracy vs Speed Tradeoff Example ===\n\n');

% Define a 2D hypercube: -1 <= x, y <= 1
A = [1 0; -1 0; 0 1; 0 -1];
b = ones(4, 1);
expected = 4.0;

fprintf('Testing 2D Hypercube (Expected Volume: 4.0)\n');
fprintf('%s\n\n', repmat('=', 1, 60));

% Test 1: Default parameters (fast, lower accuracy)
fprintf('1. Default Parameters (epsilon=1.0, walk_length=1)\n');
fprintf('   Use case: Quick estimates\n\n');
tic;
vol1 = compute_volume(A, b);
time1 = toc;
err1 = abs(vol1 - expected) / expected * 100;
fprintf('   Result: %.4f  |  Error: %.2f%%  |  Time: %.3f sec\n\n', vol1, err1, time1);

% Test 2: Medium accuracy
fprintf('2. Medium Accuracy (epsilon=0.1, walk_length=1)\n');
fprintf('   Use case: Balance between speed and accuracy\n\n');
tic;
vol2 = compute_volume(A, b, 0.1);
time2 = toc;
err2 = abs(vol2 - expected) / expected * 100;
fprintf('   Result: %.4f  |  Error: %.2f%%  |  Time: %.3f sec\n\n', vol2, err2, time2);

% Test 3: High accuracy
fprintf('3. High Accuracy (epsilon=0.01, walk_length=10)\n');
fprintf('   Use case: Research-quality results\n\n');
tic;
vol3 = compute_volume(A, b, 0.01, 10);
time3 = toc;
err3 = abs(vol3 - expected) / expected * 100;
fprintf('   Result: %.4f  |  Error: %.2f%%  |  Time: %.3f sec\n\n', vol3, err3, time3);

% Summary
fprintf('%s\n', repmat('=', 1, 60));
fprintf('Summary:\n');
fprintf('  - Lower epsilon = higher accuracy (but slower)\n');
fprintf('  - Higher walk_length = better mixing (but slower)\n');
fprintf('  - Default parameters provide good balance for exploratory analysis\n');
fprintf('  - Tighten parameters for publication-quality results\n\n');

% Recommendation
if err2 < 1.0
    fprintf('✅ RECOMMENDATION: epsilon=0.1 achieves <1%% error with reasonable speed!\n\n');
else
    fprintf('💡 TIP: Try epsilon=0.1 for better accuracy\n\n');
end
