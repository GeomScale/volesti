cmake_minimum_required(VERSION 3.15)
project(NonSymmetricConicOptimization)

set(CMAKE_CXX_STANDARD 14)

message(STATUS "CMAKE_CXX_FLAGS are ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG_CUSTOM "-O0 -fno-omit-frame-pointer  -mno-omit-leaf-frame-pointer -DIPM_USE_DOUBLE -DIPM_DOUBLE=double")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG_CUSTOM}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -DNDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -Wall -Wextra -O3 -Xclang -fopenmp -lomp -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE_AVX "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -Wall -Wextra -O3 -mavx -mfma -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE_AVX_OPENMP "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -Wall -Wextra -O3 -march=native -Xclang  -fopenmp -lomp -mavx -mfma -fno-math-errno -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double")
set(CMAKE_CXX_FLAGS_RELEASE_LONG_DOUBLE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -openmp -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=long_double")

file(GLOB TARGETS "*.h" "*.cpp" "barriers/*.h" "barriers/*.cpp")

find_package(OpenMP REQUIRED)

find_package(Python2 COMPONENTS Development NumPy)
add_executable(NonSymmetricConicOptimization ${TARGETS})

if (NOT BOOST_DIR)
    find_package(Boost COMPONENTS *boost usr/local/)
    set(BOOST_DIR ${Boost_INCLUDE_DIRS})
endif ()

if (NOT BOOST_DIR)
    message(FATAL_ERROR "This program requires the boost library, and will not be compiled. Set with flag -BOOST_DIR.")
else ()
    target_include_directories(NonSymmetricConicOptimization PRIVATE ${BOOST_DIR}
            ../../external include ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS} ../../../Eigen/eigen /usr/local/include /usr/local/opt/llvm/include)
    target_link_directories(NonSymmetricConicOptimization PRIVATE /usr/local/lib /usr/local/opt/llvm/lib /usr/local/opt/libomp/lib/)
    target_link_libraries(NonSymmetricConicOptimization Python2::Python Python2::NumPy)
endif ()