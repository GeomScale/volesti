## -*- texinfo -*-
## @deftypefn {Function File} {} test_sampling ()
## 
## Test sampling functions.
##
## @seealso{sample_points}
## @end deftypefn

function test_sampling ()
  printf ("Testing sampling functions...\n");
  
  % Test 1: Basic sampling
  printf ("  Test 1: Basic sampling...");
  P = GenCube(3);
  n = 100;
  points = sample_points(P, n);
  assert (size(points) == [3, n], "Should return correct size");
  assert (all(all(points >= -1.1)), "Points should be within bounds");
  assert (all(all(points <= 1.1)), "Points should be within bounds");
  printf (" PASSED\n");
  
  % Test 2: Different walk methods
  printf ("  Test 2: Different walk methods...");
  P = GenCube(2);
  n = 50;
  methods = {"cdhr", "rdhr", "ball", "billiard"};
  for i = 1:length(methods)
    points = sample_points(P, n, methods{i});
    assert (size(points) == [2, n], sprintf ("Method %s should work", methods{i}));
  end
  printf (" PASSED\n");
  
  % Test 3: With burn-in
  printf ("  Test 3: Sampling with burn-in...");
  P = GenCube(2);
  n = 30;
  points = sample_points(P, n, "cdhr", [], 10);
  assert (size(points) == [2, n], "Should return correct size");
  printf (" PASSED\n");
  
  % Test 4: Custom walk length
  printf ("  Test 4: Custom walk length...");
  P = GenCube(2);
  n = 20;
  points = sample_points(P, n, "cdhr", 5, 0);
  assert (size(points) == [2, n], "Should return correct size");
  printf (" PASSED\n");
  
  printf ("All sampling tests passed!\n");
endfunction

