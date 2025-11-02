## -*- texinfo -*-
## @deftypefn {Function File} {} test_volume ()
## 
## Test volume computation functions.
##
## @seealso{volume}
## @end deftypefn

function test_volume ()
  printf ("Testing volume computation...\n");
  
  % Test 1: 2D square (side length 2, should be ~4)
  printf ("  Test 1: 2D square...");
  P = GenCube(2);
  vol = volume(P);
  assert (vol > 0 && vol < 10, "Volume should be reasonable");
  printf (" PASSED (vol = %f)\n", vol);
  
  % Test 2: 3D cube (should be ~8)
  printf ("  Test 2: 3D cube...");
  P = GenCube(3);
  vol = volume(P);
  assert (vol > 0 && vol < 20, "Volume should be reasonable");
  printf (" PASSED (vol = %f)\n", vol);
  
  % Test 3: Different methods
  printf ("  Test 3: Different methods...");
  P = GenCube(2);
  vol1 = volume(P, "sequence_of_balls");
  vol2 = volume(P, "cooling_gaussians");
  assert (vol1 > 0 && vol2 > 0, "Both methods should work");
  printf (" PASSED\n");
  
  % Test 4: Custom error tolerance
  printf ("  Test 4: Custom error tolerance...");
  P = GenCube(2);
  vol = volume(P, "sequence_of_balls", 0.2);
  assert (vol > 0, "Should compute volume");
  printf (" PASSED\n");
  
  printf ("All volume tests passed!\n");
endfunction

