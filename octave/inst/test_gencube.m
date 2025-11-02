## -*- texinfo -*-
## @deftypefn {Function File} {} test_gencube ()
## 
## Test GenCube polytope generator.
##
## @seealso{GenCube}
## @end deftypefn

function test_gencube ()
  printf ("Testing GenCube function...\n");
  
  % Test 1: Default parameters
  printf ("  Test 1: Default parameters...");
  P = GenCube(3);
  assert (P.dimension == 3, "Dimension should be 3");
  assert (rows(P.A) == 6, "Should have 6 constraints");
  assert (columns(P.A) == 3, "A should have 3 columns");
  assert (all(P.b == 1), "Default scale should be 1");
  printf (" PASSED\n");
  
  % Test 2: Custom scale
  printf ("  Test 2: Custom scale...");
  P = GenCube(2, "H", 2.5);
  assert (P.dimension == 2, "Dimension should be 2");
  assert (abs(max(P.b) - 2.5) < 1e-10, "Scale should be 2.5");
  printf (" PASSED\n");
  
  % Test 3: Dimension 1
  printf ("  Test 3: Dimension 1...");
  P = GenCube(1);
  assert (P.dimension == 1, "Dimension should be 1");
  assert (rows(P.A) == 2, "Should have 2 constraints");
  printf (" PASSED\n");
  
  % Test 4: Higher dimension
  printf ("  Test 4: Higher dimension...");
  P = GenCube(5);
  assert (P.dimension == 5, "Dimension should be 5");
  assert (rows(P.A) == 10, "Should have 10 constraints");
  printf (" PASSED\n");
  
  % Test 5: Polytope structure fields
  printf ("  Test 5: Polytope structure fields...");
  P = GenCube(2);
  assert (isfield(P, "A"), "Should have field A");
  assert (isfield(P, "b"), "Should have field b");
  assert (isfield(P, "dimension"), "Should have field dimension");
  printf (" PASSED\n");
  
  printf ("All GenCube tests passed!\n");
endfunction

