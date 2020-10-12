cmake_minimum_required(VERSION 3.15)
project(NonSymmetricConicOptimization)

set(CMAKE_CXX_STANDARD 14)

set(MKLROOT /opt/intel/compilers_and_libraries_2020.1.216/mac/mkl)
set(DYLD_LIBRARY_PATH /opt/intel/compilers_and_libraries_2020.1.216/mac/tbb/lib:/opt/intel/compilers_and_libraries_2020.1.216/mac/compiler/lib:/opt/intel/compilers_and_libraries_2020.1.216/mac/mkl/lib)
set(LIBRARY_PATH /opt/intel/compilers_and_libraries_2020.1.216/mac/tbb/lib:/opt/intel/compilers_and_libraries_2020.1.216/mac/compiler/lib:/opt/intel/compilers_and_libraries_2020.1.216/mac/mkl/lib)
set(CPATH /opt/intel/compilers_and_libraries_2020.1.216/mac/mkl/include)
set(NLSPATH /opt/intel/compilers_and_libraries_2020.1.216/mac/mkl/lib/locale/%l_%t/%N)

set(CMAKE_CXX_FLAGS_DEBUG_CUSTOM "-O0 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -DIPM_USE_DOUBLE -DIPM_DOUBLE=double")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG_CUSTOM}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -DNDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE_AVX "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -Wall -Wextra -O3 -mavx -mfma -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double")

#Removed -lmkl_intel_thread and -liomp. Add again if fails.
set(MKL_OPTIONS "-DMKL_LP64 -m64 -I${MKLROOT}/include -L${DYLD_LIBRARY_PATH} -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl")

set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE_AVX_OPENMP "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -Wall -Wextra -O3 -march=native -Xclang  -fopenmp -lomp -mavx -mfma -fno-math-errno -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double -DEIGEN_USE_MKL_ALL=1 ${MKL_OPTIONS}")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE            "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -Wall -Wextra -O3 -march=native -Xclang  -fopenmp -lomp -mavx -mfma -fno-math-errno -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double ${MKL_OPTIONS}")

set(CMAKE_CXX_FLAGS_RELEASE_LONG_DOUBLE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -openmp -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=long_double")

file(GLOB TARGETS "*.h" "*.cpp" "barriers/*.h" "barriers/*.cpp"
        "../../external/ChebTools/*"
        "../../external/Padua/*"
        )

find_package(OpenMP REQUIRED)
#find_package(spdlog REQUIRED)

find_package(Python2 COMPONENTS Development NumPy)
add_executable(NonSymmetricConicOptimization ${TARGETS} utils.cpp)

if (NOT BOOST_DIR)
    find_package(Boost COMPONENTS *boost usr/local/)
    set(BOOST_DIR ${Boost_INCLUDE_DIRS})
    message("Boost found at ${Boost_INCLUDE_DIRS}")
endif ()

if (NOT BOOST_DIR)
    message(FATAL_ERROR "This program requires the boost library, and will not be compiled. Set with flag -BOOST_DIR.")
else ()
    target_include_directories(NonSymmetricConicOptimization PRIVATE ${BOOST_DIR}
            . ../../external  ${SPDLOG_DIR} include ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS}
            ../../../Eigen/eigen
            /usr/local/include
            /usr/local/opt/llvm/include)
    target_link_directories(NonSymmetricConicOptimization PRIVATE /usr/local/lib /usr/local/opt/llvm/lib /usr/local/opt/libomp/lib/ ${DYLD_LIBRARY_PATH} ${LIBRARY_PATH})
    target_link_libraries(NonSymmetricConicOptimization Python2::Python Python2::NumPy)
endif ()