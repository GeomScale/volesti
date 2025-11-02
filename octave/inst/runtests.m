## -*- texinfo -*-
## @deftypefn {Function File} {} runtests ()
## 
## Run all test suites for the VolEsti package.
##
## This function runs all individual test suites:
##   - test_gencube
##   - test_volume
##   - test_sampling
##   - volesti_test (comprehensive test suite)
##
## @seealso{test_gencube, test_volume, test_sampling, volesti_test}
## @end deftypefn

function runtests ()
  printf ("\n========================================\n");
  printf ("Running VolEsti Package Test Suites\n");
  printf ("========================================\n\n");
  
  try
    printf ("[1/4] Running GenCube tests...\n");
    test_gencube ();
    printf ("✓ GenCube tests passed\n\n");
  catch err
    printf ("✗ GenCube tests failed: %s\n\n", err.message);
  end
  
  try
    printf ("[2/4] Running volume computation tests...\n");
    test_volume ();
    printf ("✓ Volume tests passed\n\n");
  catch err
    printf ("✗ Volume tests failed: %s\n\n", err.message);
  end
  
  try
    printf ("[3/4] Running sampling tests...\n");
    test_sampling ();
    printf ("✓ Sampling tests passed\n\n");
  catch err
    printf ("✗ Sampling tests failed: %s\n\n", err.message);
  end
  
  try
    printf ("[4/4] Running comprehensive test suite...\n");
    volesti_test ();
    printf ("✓ Comprehensive tests completed\n\n");
  catch err
    printf ("✗ Comprehensive tests failed: %s\n\n", err.message);
  end
  
  printf ("========================================\n");
  printf ("All test suites completed\n");
  printf ("========================================\n\n");
endfunction

