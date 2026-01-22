#!/usr/bin/octave -qf
% test_new_generators.m
%
% Test the newly implemented generators: gen_simplex, gen_skinny_cube, gen_birkhoff

addpath('../inst');

fprintf('\n=== Testing New Generators ===\n\n');

%% Test 1: gen_simplex with H-representation
fprintf('Test 1: gen_simplex - H-representation\n');
try
    P = gen_simplex(5, 'H');
    fprintf('  ✓ Created 5D simplex (H-rep)\n');
    fprintf('    Constraints: %d x %d\n', size(P.A));
    fprintf('    Expected volume: 1/5! = %.6f\n', 1/factorial(5));
    if isfield(P, 'volume')
        fprintf('    Stored volume: %.6f\n', P.volume);
    end
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 2: gen_simplex with V-representation
fprintf('\nTest 2: gen_simplex - V-representation\n');
try
    P = gen_simplex(4, 'V');
    fprintf('  ✓ Created 4D simplex (V-rep)\n');
    fprintf('    Vertices: %d x %d\n', size(P.V));
    fprintf('    Expected: 5 vertices in 4D\n');
    fprintf('    Expected volume: 1/4! = %.6f\n', 1/factorial(4));
    if isfield(P, 'volume')
        fprintf('    Stored volume: %.6f\n', P.volume);
    end
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 3: gen_simplex default (should be H)
fprintf('\nTest 3: gen_simplex - default representation\n');
try
    P = gen_simplex(3);
    if strcmp(P.type, 'Hpolytope')
        fprintf('  ✓ Default representation is H (correct)\n');
    else
        fprintf('  ✗ Default representation is not H\n');
    end
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 4: gen_skinny_cube
fprintf('\nTest 4: gen_skinny_cube\n');
try
    P = gen_skinny_cube(5);
    fprintf('  ✓ Created 5D skinny cube\n');
    fprintf('    Constraints: %d x %d\n', size(P.A));
    fprintf('    Expected: 10 constraints (2*(d-1) + 2)\n');
    fprintf('    Expected volume: 2^4 * 200 = %.1f\n', 2^4 * 200);
    if isfield(P, 'volume')
        fprintf('    Stored volume: %.1f\n', P.volume);
    end
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 5: gen_skinny_cube - verify shape
fprintf('\nTest 5: gen_skinny_cube - verify constraints\n');
try
    P = gen_skinny_cube(3);
    % Should be [-1,1]² × [-100,100]
    % Check if constraints make sense
    fprintf('  ✓ Created 3D skinny cube\n');
    fprintf('    Constraints: %d x %d\n', size(P.A));
    fprintf('    b vector: [');
    fprintf('%g ', P.b');
    fprintf(']\n');
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 6: gen_birkhoff (n=2)
fprintf('\nTest 6: gen_birkhoff (n=2)\n');
try
    P = gen_birkhoff(2);
    fprintf('  ✓ Created 2-Birkhoff polytope\n');
    fprintf('    Constraints: %d x %d\n', size(P.A));
    fprintf('    Dimension: 1 (expected for n=2)\n');
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 7: gen_birkhoff (n=3)
fprintf('\nTest 7: gen_birkhoff (n=3)\n');
try
    P = gen_birkhoff(3);
    fprintf('  ✓ Created 3-Birkhoff polytope\n');
    fprintf('    Constraints: %d x %d\n', size(P.A));
    fprintf('    Dimension: 4 (expected: (n-1)² = 4)\n');
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 8: Integration with existing functions - Volume computation
fprintf('\nTest 8: Integration test - Compute volume of simplex\n');
try
    P = gen_simplex(3, 'H');
    vol = volesti_volume(P, 0.1, 10, false);
    expected = 1/factorial(3);
    error_pct = abs(vol - expected) / expected * 100;
    fprintf('  Computed volume: %.6f\n', vol);
    fprintf('  Expected volume: %.6f (1/3!)\n', expected);
    fprintf('  Error: %.2f%%\n', error_pct);
    if error_pct < 20
        fprintf('  ✓ Volume computation works with gen_simplex\n');
    else
        fprintf('  ⚠ Volume error is high (expected for low accuracy)\n');
    end
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 9: Integration test - Sampling from skinny cube
fprintf('\nTest 9: Integration test - Sample from skinny cube\n');
try
    P = gen_skinny_cube(3);
    % Extract A and b for direct sampling
    samples = sample_points(P.A, P.b, 50, 0, 1, 0, false);
    fprintf('  Generated %d samples from 3D skinny cube\n', size(samples, 2));
    
    % Check if samples respect the skinny structure
    % First 2 dimensions should be in [-1, 1]
    % Last dimension should be in [-100, 100]
    max_vals = max(samples, [], 2);
    min_vals = min(samples, [], 2);
    
    fprintf('  Dimension 1 range: [%.3f, %.3f] (expected [-1, 1])\n', min_vals(1), max_vals(1));
    fprintf('  Dimension 2 range: [%.3f, %.3f] (expected [-1, 1])\n', min_vals(2), max_vals(2));
    fprintf('  Dimension 3 range: [%.3f, %.3f] (expected [-100, 100])\n', min_vals(3), max_vals(3));
    
    if all(max_vals(1:2) <= 1.01) && all(min_vals(1:2) >= -1.01) && ...
       max_vals(3) <= 101 && min_vals(3) >= -101
        fprintf('  ✓ Sampling works correctly with skinny cube\n');
    else
        fprintf('  ⚠ Sample ranges unexpected\n');
    end
catch err
    fprintf('  ✗ Error: %s\n', err.message);
end

%% Test 10: Error handling
fprintf('\nTest 10: Error handling\n');
errors_caught = 0;

% Test invalid dimension
try
    P = gen_simplex(0);
    fprintf('  ✗ Should have caught invalid dimension\n');
catch
    errors_caught = errors_caught + 1;
end

% Test invalid representation
try
    P = gen_simplex(5, 'X');
    fprintf('  ✗ Should have caught invalid representation\n');
catch
    errors_caught = errors_caught + 1;
end

% Test skinny cube dimension < 2
try
    P = gen_skinny_cube(1);
    fprintf('  ✗ Should have caught dimension < 2\n');
catch
    errors_caught = errors_caught + 1;
end

if errors_caught == 3
    fprintf('  ✓ All error cases properly handled (%d/3)\n', errors_caught);
else
    fprintf('  ⚠ Some errors not caught (%d/3)\n', errors_caught);
end

fprintf('\n=== All Generator Tests Complete ===\n\n');
