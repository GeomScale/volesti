cmake_minimum_required(VERSION 3.15)
project(NonSymmetricConicOptimization)

set(CMAKE_CXX_STANDARD 14)
set(TARGETS NonSymmetricIPM.cpp  NonSymmetricIPM.h LHSCB.cpp LHSCB.h MonomialsClass.cpp MonomialsClass.h
        EnvelopeProblemSOS.cpp EnvelopeProblemSDP.cpp EnvelopeProblemSDP.h utils.cpp EnvelopeProblemSOS.h)

find_package(Python2 COMPONENTS Development NumPy)

add_executable(NonSymmetricConicOptimization ${TARGETS} sos.cpp)
target_include_directories(NonSymmetricConicOptimization PRIVATE ../../external include ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS})
target_link_libraries(NonSymmetricConicOptimization Python2::Python Python2::NumPy)
