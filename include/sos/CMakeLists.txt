cmake_minimum_required(VERSION 3.15)
project(NonSymmetricConicOptimization)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "-DIPM_USE_DOUBLE")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -mavx -DNDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE  "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -mavx -DNDEBUG=1 -DIPM_USE_DOUBLE")

set(TARGETS NonSymmetricIPM.cpp  NonSymmetricIPM.h LHSCB.cpp LHSCB.h MonomialsClass.cpp MonomialsClass.h
        EnvelopeProblemSOS.cpp EnvelopeProblemSDP.cpp EnvelopeProblemSDP.h utils.cpp EnvelopeProblemSOS.h)

find_package(Python2 COMPONENTS Development NumPy)

link_directories(${CMAKE_SOURCE_DIR}/lib /usr/local/lib)
add_executable(NonSymmetricConicOptimization ${TARGETS} sos.cpp)

# TODO: Add used Boost Header files into external folder
target_include_directories(NonSymmetricConicOptimization PRIVATE /usr/local/Cellar/libomp/10.0.0/include /usr/local/Cellar/boost/1.72.0/include ../../external include ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS}
        )

target_link_libraries(NonSymmetricConicOptimization Python2::Python Python2::NumPy )
