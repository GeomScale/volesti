## -*- texinfo -*-
## @deftypefn {Script File} {} volesti_example ()
## 
## Run example demonstrating VolEsti functionality.
##
## This script demonstrates:
##   - Generating a hypercube polytope
##   - Computing its volume
##   - Sampling uniform points from it
##
## @seealso{volume, sample_points, GenCube}
## @end deftypefn

function volesti_example ()
  printf ("\n=== VolEsti Octave Package Example ===\n\n");
  
  % Generate a 3-dimensional cube
  printf ("Generating a 3-dimensional cube...\n");
  P = GenCube(3);
  printf ("Done. Dimension: %d\n", P.dimension);
  printf ("Number of constraints: %d\n\n", rows(P.A));
  
  % Compute volume
  printf ("Computing volume (this may take a moment)...\n");
  tic;
  vol = volume(P);
  elapsed = toc;
  printf ("Volume: %f (exact: 8.0)\n", vol);
  printf ("Computation time: %.2f seconds\n\n", elapsed);
  
  % Sample points
  n_points = 100;
  printf ("Sampling %d uniform points...\n", n_points);
  tic;
  points = sample_points(P, n_points);
  elapsed = toc;
  printf ("Sampling time: %.2f seconds\n", elapsed);
  printf ("Sampled points shape: %d x %d\n\n", rows(points), columns(points));
  
  % Plot if dimension allows
  if (P.dimension >= 2)
    printf ("Plotting first two dimensions of sampled points...\n");
    figure;
    plot(points(1,:), points(2,:), 'o');
    xlabel('x_1');
    ylabel('x_2');
    title('Uniform samples from 3D cube (projected to first two dimensions)');
    grid on;
    printf ("Plot displayed.\n\n");
  endif
  
  printf ("Example completed successfully!\n\n");
endfunction

