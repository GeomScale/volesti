# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build

# Include any dependencies generated for this target.
include CMakeFiles/boundary_oracles_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/boundary_oracles_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/boundary_oracles_test.dir/flags.make

CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.o: CMakeFiles/boundary_oracles_test.dir/flags.make
CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.o: ../boundary_oracles_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.o -c /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/boundary_oracles_test.cpp

CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/boundary_oracles_test.cpp > CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.i

CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/boundary_oracles_test.cpp -o CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.s

# Object files for target boundary_oracles_test
boundary_oracles_test_OBJECTS = \
"CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.o"

# External object files for target boundary_oracles_test
boundary_oracles_test_EXTERNAL_OBJECTS = \
"/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles/test_main.dir/test_main.cpp.o"

boundary_oracles_test: CMakeFiles/boundary_oracles_test.dir/boundary_oracles_test.cpp.o
boundary_oracles_test: CMakeFiles/test_main.dir/test_main.cpp.o
boundary_oracles_test: CMakeFiles/boundary_oracles_test.dir/build.make
boundary_oracles_test: /usr/lib/lp_solve/liblpsolve55.so
boundary_oracles_test: CMakeFiles/boundary_oracles_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable boundary_oracles_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/boundary_oracles_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/boundary_oracles_test.dir/build: boundary_oracles_test

.PHONY : CMakeFiles/boundary_oracles_test.dir/build

CMakeFiles/boundary_oracles_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/boundary_oracles_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/boundary_oracles_test.dir/clean

CMakeFiles/boundary_oracles_test.dir/depend:
	cd /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles/boundary_oracles_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/boundary_oracles_test.dir/depend
