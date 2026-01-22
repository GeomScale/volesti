# Makefile for Volesti Octave Interface
# Compiles the C++ wrappers into .oct plugins using mkoctfile

# Compiler
MKOCTFILE = mkoctfile

# Volesti include path (parent directory of include/)
VOLESTI_ROOT = ../..
VOLESTI_INCLUDE = $(VOLESTI_ROOT)/include

# External dependencies (Eigen, Boost, Spectra, lpsolve) - from CMake FetchContent
# IMPORTANT: Before running make, you must run CMake from the repository root
# (e.g., cmake -S . -B build) so FetchContent downloads and populates the
# external/_deps tree with Eigen, Boost, and lpsolve dependencies.
EXTERNAL = $(VOLESTI_ROOT)/external
EIGEN_INCLUDE = $(EXTERNAL)/_deps/eigen-src
BOOST_INCLUDE = $(EXTERNAL)/_deps/boost-src
LPSOLVE_INCLUDE = $(EXTERNAL)/_deps/lpsolve-src

# Compiler flags
# Note: EXTERNAL is added so includes like "Spectra/include/..." work
CXXFLAGS = -I$(VOLESTI_INCLUDE) -I$(EXTERNAL) -I$(EIGEN_INCLUDE) -I$(BOOST_INCLUDE) -I$(LPSOLVE_INCLUDE) -std=c++17 -O2

# Linker flags - use shared library and set rpath for runtime loading
LDFLAGS = -L/usr/lib/lp_solve -llpsolve55 -Wl,-rpath,/usr/lib/lp_solve

# Default target: build all .oct plugins
all: inst/volume.oct inst/sample_points.oct

# Build the volume.oct plugin
inst/volume.oct: src/volume.cpp
	$(MKOCTFILE) $(CXXFLAGS) src/volume.cpp -o inst/volume.oct $(LDFLAGS)

# Build the sample_points.oct plugin
inst/sample_points.oct: src/sample_points.cpp
	$(MKOCTFILE) $(CXXFLAGS) src/sample_points.cpp -o inst/sample_points.oct $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f inst/volume.oct inst/sample_points.oct src/*.o

# Run all tests (like ctest)
test: inst/volume.oct inst/sample_points.oct
	@echo "Running Octave binding tests..."
	@octave --eval "addpath('inst'); cd('test'); test_Hvol" && \
	octave --eval "addpath('inst'); cd('test'); test_Vvol" && \
	octave --eval "addpath('inst'); cd('test'); test_generators" && \
	echo "All tests passed!" || echo "Some tests failed!"

# Phony targets
.PHONY: all clean test
