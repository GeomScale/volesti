@echo off
REM Build script for VolEsti Octave MEX files (Windows)

REM Create build directory
if not exist build mkdir build
cd build

REM Configure with CMake
cmake ..

REM Build
cmake --build . --config Release

REM Check if build succeeded
if %ERRORLEVEL% EQU 0 (
    echo Build successful! MEX files are in ..\inst\
    echo To install, run from Octave: pkg install ..\
) else (
    echo Build failed!
    exit /b 1
)

cd ..

