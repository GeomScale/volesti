# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/doyouseeme/volesti/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/doyouseeme/volesti/test

# Include any dependencies generated for this target.
include CMakeFiles/volume_cg_hpolytope.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/volume_cg_hpolytope.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/volume_cg_hpolytope.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/volume_cg_hpolytope.dir/flags.make

CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o: CMakeFiles/volume_cg_hpolytope.dir/flags.make
CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o: volume_cg_hpolytope.cpp
CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o: CMakeFiles/volume_cg_hpolytope.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/doyouseeme/volesti/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o -MF CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o.d -o CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o -c /home/doyouseeme/volesti/test/volume_cg_hpolytope.cpp

CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/doyouseeme/volesti/test/volume_cg_hpolytope.cpp > CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.i

CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/doyouseeme/volesti/test/volume_cg_hpolytope.cpp -o CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.s

# Object files for target volume_cg_hpolytope
volume_cg_hpolytope_OBJECTS = \
"CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o"

# External object files for target volume_cg_hpolytope
volume_cg_hpolytope_EXTERNAL_OBJECTS = \
"/home/doyouseeme/volesti/test/CMakeFiles/test_main.dir/test_main.cpp.o"

volume_cg_hpolytope: CMakeFiles/volume_cg_hpolytope.dir/volume_cg_hpolytope.cpp.o
volume_cg_hpolytope: CMakeFiles/test_main.dir/test_main.cpp.o
volume_cg_hpolytope: CMakeFiles/volume_cg_hpolytope.dir/build.make
volume_cg_hpolytope: liblp_solve.a
volume_cg_hpolytope: CMakeFiles/volume_cg_hpolytope.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/doyouseeme/volesti/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable volume_cg_hpolytope"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/volume_cg_hpolytope.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/volume_cg_hpolytope.dir/build: volume_cg_hpolytope
.PHONY : CMakeFiles/volume_cg_hpolytope.dir/build

CMakeFiles/volume_cg_hpolytope.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/volume_cg_hpolytope.dir/cmake_clean.cmake
.PHONY : CMakeFiles/volume_cg_hpolytope.dir/clean

CMakeFiles/volume_cg_hpolytope.dir/depend:
	cd /home/doyouseeme/volesti/test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/doyouseeme/volesti/test /home/doyouseeme/volesti/test /home/doyouseeme/volesti/test /home/doyouseeme/volesti/test /home/doyouseeme/volesti/test/CMakeFiles/volume_cg_hpolytope.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/volume_cg_hpolytope.dir/depend

