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
include CMakeFiles/volume_sob_hpolytope.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/volume_sob_hpolytope.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/volume_sob_hpolytope.dir/flags.make

CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.o: CMakeFiles/volume_sob_hpolytope.dir/flags.make
CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.o: ../volume_sob_hpolytope.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.o -c /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/volume_sob_hpolytope.cpp

CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/volume_sob_hpolytope.cpp > CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.i

CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/volume_sob_hpolytope.cpp -o CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.s

# Object files for target volume_sob_hpolytope
volume_sob_hpolytope_OBJECTS = \
"CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.o"

# External object files for target volume_sob_hpolytope
volume_sob_hpolytope_EXTERNAL_OBJECTS = \
"/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles/test_main.dir/test_main.cpp.o"

volume_sob_hpolytope: CMakeFiles/volume_sob_hpolytope.dir/volume_sob_hpolytope.cpp.o
volume_sob_hpolytope: CMakeFiles/test_main.dir/test_main.cpp.o
volume_sob_hpolytope: CMakeFiles/volume_sob_hpolytope.dir/build.make
volume_sob_hpolytope: /usr/lib/lp_solve/liblpsolve55.so
volume_sob_hpolytope: CMakeFiles/volume_sob_hpolytope.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable volume_sob_hpolytope"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/volume_sob_hpolytope.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/volume_sob_hpolytope.dir/build: volume_sob_hpolytope

.PHONY : CMakeFiles/volume_sob_hpolytope.dir/build

CMakeFiles/volume_sob_hpolytope.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/volume_sob_hpolytope.dir/cmake_clean.cmake
.PHONY : CMakeFiles/volume_sob_hpolytope.dir/clean

CMakeFiles/volume_sob_hpolytope.dir/depend:
	cd /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build /home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/CMakeFiles/volume_sob_hpolytope.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/volume_sob_hpolytope.dir/depend

