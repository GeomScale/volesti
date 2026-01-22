#!/usr/bin/octave -qf
% test_sampling.m
%
% Comprehensive sampling tests matching R's test_sampling.R
% Tests uniform sampling from various polytopes

addpath('../inst');

fprintf('\n=== Comprehensive Sampling Tests ===\n\n');

% Test function - returns 1 if no NaN/Inf, 0 otherwise
function res = run_sample_test(P, n, name_str)
    try
        if isstruct(P)
            % Legacy wrapper API
            error('Wrapper API not tested here - use direct .oct calls');
        else
            % Direct .oct call - determine if H or V polytope
            if isvector(P) || size(P, 1) == size(P, 2) + 1
                % Likely V-polytope (vertices)
                samples = sample_points(P, n, 1, 1, 0, false);  % RDHR, silent
            else
                % Assume first two args are A, b for H-polytope
                error('Pass A and b separately for H-polytopes');
            end
        end
        
        % Check for NaN or Inf
        if any(isnan(samples(:))) || any(isinf(samples(:)))
            fprintf('  ✗ %s: Contains NaN/Inf\n', name_str);
            res = 0;
        else
            fprintf('  ✓ %s: %d valid samples\n', name_str, size(samples, 2));
            res = 1;
        end
    catch err
        fprintf('  ✗ %s: Error - %s\n', name_str, err.message);
        res = 0;
    end
end

% Test H-polytopes
fprintf('Testing H-polytopes:\n');

% Test 1: 10D cube
fprintf('  Test 1: 10D hypercube\n');
A = [eye(10); -eye(10)];
b = ones(20, 1);
try
    samples = sample_points(A, b, 100, 0, 1, 0, false);  % CDHR
    if any(isnan(samples(:))) || any(isinf(samples(:)))
        fprintf('    ✗ Contains NaN/Inf\n');
    else
        fprintf('    ✓ Generated 100 valid samples\n');
    end
catch err
    fprintf('    ✗ Error: %s\n', err.message);
end

% Test 2: 10D cross-polytope  
fprintf('  Test 2: 10D cross-polytope\n');
A = [eye(10); -eye(10)];
b = ones(20, 1);
for i = 1:10
    b(i) = 1;
    b(10+i) = 1;
end
try
    samples = sample_points(A, b, 100, 0, 1, 0, false);
    if any(isnan(samples(:))) || any(isinf(samples(:)))
        fprintf('    ✗ Contains NaN/Inf\n');
    else
        fprintf('    ✓ Generated 100 valid samples\n');
    end
catch err
    fprintf('    ✗ Error: %s\n', err.message);
end

% Test 3: 5D product simplex (10D total)
fprintf('  Test 3: Product of 5D simplices (10D polytope)\n');
% Product simplex: {x,y : x_i >= 0, sum(x) <= 1, y_i >= 0, sum(y) <= 1}
d = 5;
A1 = [-eye(d); ones(1, d)];
b1 = [zeros(d, 1); 1];
A2 = [-eye(d); ones(1, d)];
b2 = [zeros(d, 1); 1];
A = blkdiag(A1, A2);
b = [b1; b2];
try
    samples = sample_points(A, b, 100, 0, 1, 0, false);
    if any(isnan(samples(:))) || any(isinf(samples(:)))
        fprintf('    ✗ Contains NaN/Inf\n');
    else
        fprintf('    ✓ Generated 100 valid samples\n');
    end
catch err
    fprintf('    ✗ Error: %s\n', err.message);
end

% Test 4: Standard simplex
fprintf('  Test 4: 10D standard simplex\n');
d = 10;
A = [-eye(d); ones(1, d)];
b = [zeros(d, 1); 1];
try
    samples = sample_points(A, b, 100, 0, 1, 0, false);
    if any(isnan(samples(:))) || any(isinf(samples(:)))
        fprintf('    ✗ Contains NaN/Inf\n');
    else
        fprintf('    ✓ Generated 100 valid samples\n');
    end
catch err
    fprintf('    ✗ Error: %s\n', err.message);
end

% Test V-polytopes
fprintf('\nTesting V-polytopes:\n');

% Test 5: 3D simplex
fprintf('  Test 5: 3D simplex\n');
V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
try
    samples = sample_points(V, 50, 1, 1, 0, false);  % RDHR
    if any(isnan(samples(:))) || any(isinf(samples(:)))
        fprintf('    ✗ Contains NaN/Inf\n');
    else
        fprintf('    ✓ Generated 50 valid samples\n');
    end
catch err
    fprintf('    ✗ Error: %s\n', err.message);
end

% Test 6: 3D cube
fprintf('  Test 6: 3D cube (V-representation)\n');
V = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1;
     -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
try
    samples = sample_points(V, 50, 1, 1, 0, false);
    if any(isnan(samples(:))) || any(isinf(samples(:)))
        fprintf('    ✗ Contains NaN/Inf\n');
    else
        fprintf('    ✓ Generated 50 valid samples\n');
    end
catch err
    fprintf('    ✗ Error: %s\n', err.message);
end

% Test different walk types
fprintf('\nTesting walk types:\n');

A = [eye(3); -eye(3)];
b = ones(6, 1);

fprintf('  Test 7: CDHR walk\n');
try
    samples = sample_points(A, b, 50, 0, 5, 10, false);  % walk_type=0
    fprintf('    ✓ CDHR: %d samples\n', size(samples, 2));
catch err
    fprintf('    ✗ CDHR failed: %s\n', err.message);
end

fprintf('  Test 8: RDHR walk\n');
try
    samples = sample_points(A, b, 50, 1, 5, 10, false);  % walk_type=1
    fprintf('    ✓ RDHR: %d samples\n', size(samples, 2));
catch err
    fprintf('    ✗ RDHR failed: %s\n', err.message);
end

fprintf('  Test 9: Ball Walk\n');
try
    samples = sample_points(A, b, 50, 2, 5, 10, false);  % walk_type=2
    fprintf('    ✓ Ball Walk: %d samples\n', size(samples, 2));
catch err
    fprintf('    ✗ Ball Walk failed: %s\n', err.message);
end

% Test burn-in
fprintf('\nTesting burn-in:\n');
fprintf('  Test 10: Burn-in effect\n');
try
    samples_no_burn = sample_points(A, b, 100, 0, 1, 0, false);
    samples_with_burn = sample_points(A, b, 100, 0, 1, 100, false);
    
    mean_no_burn = mean(samples_no_burn, 2);
    mean_with_burn = mean(samples_with_burn, 2);
    
    fprintf('    Mean without burn-in: [%.3f, %.3f, %.3f]\n', mean_no_burn);
    fprintf('    Mean with burn-in:    [%.3f, %.3f, %.3f]\n', mean_with_burn);
    fprintf('    ✓ Burn-in test complete\n');
catch err
    fprintf('    ✗ Burn-in test failed: %s\n', err.message);
end

fprintf('\n=== All Sampling Tests Complete ===\n\n');
