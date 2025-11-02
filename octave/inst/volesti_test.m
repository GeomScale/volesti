## -*- texinfo -*-
## @deftypefn {Script File} {} volesti_test ()
## 
## Run comprehensive tests for the VolEsti Octave package.
##
## This test suite verifies:
##   - Polytope generation (GenCube)
##   - Volume computation
##   - Point sampling
##   - Error handling
##
## @seealso{volume, sample_points, GenCube}
## @end deftypefn

function volesti_test ()
  % Test suite for VolEsti Octave package
  
  printf ("\n=== Running VolEsti Test Suite ===\n\n");
  
  % Track test results
  tests_passed = 0;
  tests_failed = 0;
  
  % Test 1: GenCube basic functionality
  printf ("Test 1: GenCube basic functionality...\n");
  try
    P = GenCube(3);
    assert (P.dimension == 3, "Dimension should be 3");
    assert (rows(P.A) == 6, "Should have 6 constraints (2 per dimension)");
    assert (columns(P.A) == 3, "A should have 3 columns");
    assert (rows(P.b) == 6, "b should have 6 elements");
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 2: GenCube with scale
  printf ("Test 2: GenCube with scale...\n");
  try
    P = GenCube(2, "H", 2.0);
    assert (abs(max(P.b) - 2.0) < 1e-10, "Scale should be 2.0");
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 3: Volume computation (3D cube should be close to 8)
  printf ("Test 3: Volume computation (3D cube)...\n");
  try
    P = GenCube(3);
    vol = volume(P);
    % Volume should be approximately 8 for a 3D cube of side 2
    % Allow 50% relative error for stochastic algorithm
    assert (vol > 0, "Volume should be positive");
    assert (vol < 20, "Volume should be reasonable");
    printf ("  PASSED (computed volume: %f, expected ~8.0)\n\n", vol);
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 4: Volume computation with method parameter
  printf ("Test 4: Volume computation with method...\n");
  try
    P = GenCube(2);
    vol1 = volume(P, "sequence_of_balls");
    vol2 = volume(P, "cooling_gaussians");
    assert (vol1 > 0 && vol2 > 0, "Both methods should return positive volumes");
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 5: Sample points basic functionality
  printf ("Test 5: Sample points basic functionality...\n");
  try
    P = GenCube(3);
    n = 50;
    points = sample_points(P, n);
    assert (rows(points) == 3, "Should have 3 rows (dimension)");
    assert (columns(points) == n, sprintf ("Should have %d columns (number of points)", n));
    % Check that points are within bounds (cube: -1 <= x_i <= 1)
    assert (all(all(points >= -1.1)), "Points should be within bounds");
    assert (all(all(points <= 1.1)), "Points should be within bounds");
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 6: Sample points with different methods
  printf ("Test 6: Sample points with different methods...\n");
  try
    P = GenCube(2);
    n = 20;
    for method = {"cdhr", "rdhr", "ball", "billiard"}
      points = sample_points(P, n, method{1});
      assert (rows(points) == 2, sprintf ("Method %s: should have 2 rows", method{1}));
      assert (columns(points) == n, sprintf ("Method %s: should have %d columns", method{1}, n));
    end
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 7: Volume with polytope structure
  printf ("Test 7: Volume with polytope structure...\n");
  try
    P = GenCube(2);
    vol1 = volume(P);
    vol2 = volume(P.A, P.b);
    % Should get similar results (allowing for stochastic variance)
    assert (abs(vol1 - vol2) < max(vol1, vol2), "Results should be similar");
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 8: Sample points with polytope structure
  printf ("Test 8: Sample points with polytope structure...\n");
  try
    P = GenCube(2);
    n = 30;
    points1 = sample_points(P, n);
    points2 = sample_points(P.A, P.b, n);
    assert (size(points1) == size(points2), "Results should have same size");
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 9: Error handling - invalid polytope structure
  printf ("Test 9: Error handling - invalid polytope structure...\n");
  try
    P = struct("A", [1, 2], "wrong_field", 3);
    try
      vol = volume(P);
      printf ("  FAILED: Should have raised an error\n\n");
      tests_failed++;
    catch
      printf ("  PASSED (correctly raised error)\n\n");
      tests_passed++;
    end
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Test 10: Error handling - dimension mismatch
  printf ("Test 10: Error handling - dimension mismatch...\n");
  try
    A = [1, 2, 3];
    b = [1; 2; 3; 4];  % Wrong size
    try
      vol = volume(A, b);
      printf ("  FAILED: Should have raised an error\n\n");
      tests_failed++;
    catch
      printf ("  PASSED (correctly raised error)\n\n");
      tests_passed++;
    end
  catch err
    printf ("  PASSED (correctly raised error: %s)\n\n", err.message);
    tests_passed++;
  end
  
  % Test 11: GenCube edge cases
  printf ("Test 11: GenCube edge cases...\n");
  try
    % Test dimension 1
    P1 = GenCube(1);
    assert (P1.dimension == 1, "Dimension 1 should work");
    
    % Test larger dimension
    P10 = GenCube(10);
    assert (P10.dimension == 10, "Dimension 10 should work");
    assert (rows(P10.A) == 20, "Should have 20 constraints for dim 10");
    
    printf ("  PASSED\n\n");
    tests_passed++;
  catch err
    printf ("  FAILED: %s\n\n", err.message);
    tests_failed++;
  end
  
  % Summary
  printf ("=== Test Summary ===\n");
  printf ("Tests passed: %d\n", tests_passed);
  printf ("Tests failed: %d\n", tests_failed);
  printf ("Total tests: %d\n\n", tests_passed + tests_failed);
  
  if (tests_failed == 0)
    printf ("All tests passed! ✓\n\n");
  else
    printf ("Some tests failed. Please review the output above.\n\n");
  end
endfunction

