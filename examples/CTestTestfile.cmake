# CMake generated Testfile for 
# Source directory: /home/tolis/Dropbox/R-project/volume_approximation/examples
# Build directory: /home/tolis/Dropbox/R-project/volume_approximation/examples
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(compilation_of__vol "/usr/bin/cmake" "--build" "/home/tolis/Dropbox/R-project/volume_approximation/examples" "--target" "vol")
set_tests_properties(compilation_of__vol PROPERTIES  FIXTURES_SETUP "vol" LABELS "GeomRandWalks_example")
add_test(execution___of__vol "/home/tolis/Dropbox/R-project/volume_approximation/examples/vol")
set_tests_properties(execution___of__vol PROPERTIES  DEPENDS "compilation_of__vol" FIXTURES_REQUIRED "GeomRandWalks_example;vol" LABELS "GeomRandWalks_example" WORKING_DIRECTORY "/home/tolis/Dropbox/R-project/volume_approximation/examples/__exec_test_dir")
add_test(GeomRandWalks_example_SetupFixture "/usr/bin/cmake" "-E" "copy_directory" "/home/tolis/Dropbox/R-project/volume_approximation/examples" "/home/tolis/Dropbox/R-project/volume_approximation/examples/__exec_test_dir")
set_tests_properties(GeomRandWalks_example_SetupFixture PROPERTIES  FIXTURES_SETUP "GeomRandWalks_example" LABELS "GeomRandWalks_example")
add_test(GeomRandWalks_example_CleanupFixture "/usr/bin/cmake" "-E" "remove_directory" "/home/tolis/Dropbox/R-project/volume_approximation/examples/__exec_test_dir")
set_tests_properties(GeomRandWalks_example_CleanupFixture PROPERTIES  FIXTURES_CLEANUP "GeomRandWalks_example" LABELS "GeomRandWalks_example")
