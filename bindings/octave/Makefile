# Makefile for Volesti Octave Interface
# Compiles the C++ wrapper into a .oct plugin using mkoctfile

# Compiler
MKOCTFILE = mkoctfile

# Volesti include path (parent directory of include/)
VOLESTI_ROOT = ../..
VOLESTI_INCLUDE = $(VOLESTI_ROOT)/include

# External dependencies (Eigen, Boost, Spectra) - from CMake FetchContent
# IMPORTANT: Before running make, you must run CMake from the repository root
# (e.g., cmake -S . -B build) so FetchContent downloads and populates the
# external/_deps tree with Eigen and Boost dependencies.
EXTERNAL = $(VOLESTI_ROOT)/external
EIGEN_INCLUDE = $(EXTERNAL)/_deps/eigen-src
BOOST_INCLUDE = $(EXTERNAL)/_deps/boost-src

# Compiler flags
# Note: EXTERNAL is added so includes like "Spectra/include/..." work
# DISABLE_LPSOLVE: PoC doesn't need LP solver functionality
CXXFLAGS = -I$(VOLESTI_INCLUDE) -I$(EXTERNAL) -I$(EIGEN_INCLUDE) -I$(BOOST_INCLUDE) -std=c++17 -O2 -DDISABLE_LPSOLVE

# Source and output
SRC = src/volume.cpp
TARGET = inst/volume.oct

# Default target: build the .oct plugin
all: $(TARGET)

# Build the .oct plugin from the C++ wrapper
$(TARGET): $(SRC)
	$(MKOCTFILE) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Clean build artifacts
clean:
	rm -f $(TARGET) src/*.o

# Phony targets
.PHONY: all clean
