#!/bin/bash
# Build script for VolEsti Octave MEX files

# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake ..

# Build
make

# Check if build succeeded
if [ $? -eq 0 ]; then
    echo "Build successful! MEX files are in ../inst/"
    echo "To install, run from Octave: pkg install ../"
else
    echo "Build failed!"
    exit 1
fi

