cmake_minimum_required(VERSION 3.15)
project(NonSymmetricConicOptimization)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "-DIPM_USE_DOUBLE")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -DNDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -openmp -DNDEBUG=1 -DIPM_USE_DOUBLE")
set(CMAKE_CXX_FLAGS_RELEASE_LONG_DOUBLE  "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -openmp -DNDEBUG=1 -DIPM_USE_DOUBLE -DIPM_DOUBLE=long_double")



set(TARGETS NonSymmetricIPM.cpp NonSymmetricIPM.h LHSCB.cpp LHSCB.h MonomialsClass.cpp MonomialsClass.h
        EnvelopeProblemSOS.cpp EnvelopeProblemSDP.cpp EnvelopeProblemSDP.h utils.cpp EnvelopeProblemSOS.h)

find_package(Python2 COMPONENTS Development NumPy)
add_executable(NonSymmetricConicOptimization ${TARGETS} sos.cpp)

if (NOT BOOST_DIR)
    find_package(Boost COMPONENTS *boost usr/local/Cellar/boost/1.72.0/include)
    set(BOOST_DIR ${Boost_INCLUDE_DIRS})
endif()

if (NOT BOOST_DIR)
    message(FATAL_ERROR "This program requires the boost library, and will not be compiled. Set with flag -BOOST_DIR.")
else ()
    target_include_directories(NonSymmetricConicOptimization PRIVATE ${BOOST_DIR}
            ../../external include ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS})

    target_link_libraries(NonSymmetricConicOptimization Python2::Python Python2::NumPy)
endif()