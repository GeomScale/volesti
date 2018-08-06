# CMake generated Testfile for 
# Source directory: /home/tolis/volume_approximation/test
# Build directory: /home/tolis/volume_approximation/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(volume_cube "/home/tolis/volume_approximation/test/volume_test" "-tc=cube")
add_test(volume_cross "/home/tolis/volume_approximation/test/volume_test" "-tc=cross")
add_test(volume_birkhoff "/home/tolis/volume_approximation/test/volume_test" "-tc=birkhoff")
add_test(volume_prod_simplex "/home/tolis/volume_approximation/test/volume_test" "-tc=prod_simplex")
add_test(volume_skinny_cube "/home/tolis/volume_approximation/test/volume_test" "-tc=skinny_cube")
add_test(volumeCV_cube "/home/tolis/volume_approximation/test/volumeCV_test" "-tc=cube")
add_test(volumeCV_cross "/home/tolis/volume_approximation/test/volumeCV_test" "-tc=cross")
add_test(volumeCV_birkhoff "/home/tolis/volume_approximation/test/volumeCV_test" "-tc=birkhoff")
add_test(volumeCV_prod_simplex "/home/tolis/volume_approximation/test/volumeCV_test" "-tc=prod_simplex")
add_test(volume_float_cube "/home/tolis/volume_approximation/test/float_VE" "-tc=cube")
add_test(volume_float_cross "/home/tolis/volume_approximation/test/float_VE" "-tc=cross")
add_test(volume_float_birkhoff "/home/tolis/volume_approximation/test/float_VE" "-tc=birkhoff")
add_test(volume_float_prod_simplex "/home/tolis/volume_approximation/test/float_VE" "-tc=prod_simplex")
add_test(volume_float_skinny_cube "/home/tolis/volume_approximation/test/float_VE" "-tc=skinny_cube")
add_test(volume_long_double_cube "/home/tolis/volume_approximation/test/long_double_CV" "-tc=cube")
add_test(volumeCV_long_double_cross "/home/tolis/volume_approximation/test/long_double_CV" "-tc=cross")
add_test(volumeCV_long_double_birkhoff "/home/tolis/volume_approximation/test/long_double_CV" "-tc=birkhoff")
add_test(volumeCV_long_double_prod_simplex "/home/tolis/volume_approximation/test/long_double_CV" "-tc=prod_simplex")
add_test(cheb_cube "/home/tolis/volume_approximation/test/cheb_test" "-tc=cheb_cube")
add_test(cheb_cross "/home/tolis/volume_approximation/test/cheb_test" "-tc=cheb_cross")
add_test(cheb_birkhoff "/home/tolis/volume_approximation/test/cheb_test" "-tc=cheb_birkhoff")
add_test(cheb_prod_simplex "/home/tolis/volume_approximation/test/cheb_test" "-tc=cheb_prod_simplex")
add_test(cheb_simplex "/home/tolis/volume_approximation/test/cheb_test" "-tc=cheb_simplex")
add_test(cheb_skinny_cube "/home/tolis/volume_approximation/test/cheb_test" "-tc=cheb_skinny_cube")
add_test(round_skinny_cube "/home/tolis/volume_approximation/test/rounding_test" "-tc=round_skinny_cube")
add_test(round_rot_skinny_cube "/home/tolis/volume_approximation/test/rounding_test" "-tc=round_rot_skinny_cube")
