# Makefile for Volesti Octave Interface
# Compiles the C++ wrapper into a .oct plugin using mkoctfile

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

# Source and output
SRC = src/volume.cpp
TARGET = inst/volume.oct

# Linker flags - use shared library and set rpath for runtime loading
LDFLAGS = -L/usr/lib/lp_solve -llpsolve55 -Wl,-rpath,/usr/lib/lp_solve

# Default target: build the .oct plugin
all: $(TARGET)

# Build the .oct plugin from the C++ wrapper
$(TARGET): $(SRC)
	$(MKOCTFILE) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(TARGET) src/*.o

# Run all tests (like ctest)
test: $(TARGET)
	@echo "Running Octave binding tests..."
	@octave --eval "addpath('inst'); cd('test'); test_Hvol" && \
	octave --eval "addpath('inst'); cd('test'); test_Vvol" && \
	echo "All tests passed!" || echo "Some tests failed!"

# Phony targets
.PHONY: all clean test
