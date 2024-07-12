cmake_minimum_required(VERSION 3.15)
project(EnvelopeProblem)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_DEBUG_CUSTOM "-O0 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -DIPM_USE_DOUBLE -DIPM_DOUBLE=double")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG_CUSTOM}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -DNDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE_AVX "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -Wall -Wextra -O3 -mavx -mfma -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double")

#Environment variables MKLROOT, DYLD_LIBRARY_PATH, LIBRARY_PATH, CPATH and NLSPATH must be set to use MKL. See more information here: https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
set(MKL_OPTIONS "-DMKL_LP64 -m64 -I$ENV{MKLROOT}/include")

set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE_AVX_OPENMP "${CMAKE_CXX_FLAGS_DEBUG} -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -O3 -march=native -Xclang  -fopenmp -mavx -mfma -fno-math-errno -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double -DEIGEN_USE_MKL_ALL=1 ${MKL_OPTIONS}")
set(CMAKE_CXX_FLAGS_RELEASE_DOUBLE "${CMAKE_CXX_FLAGS} -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -Wall -Wextra -O3 -march=native -Xclang  -fopenmp -mavx -mfma -fno-math-errno -DNDEBUG=1 -DIPM_USE_DOUBLE -DDIGITS_PRECISION=250 -DIPM_DOUBLE=double ${MKL_OPTIONS}")

file(GLOB TARGETS "*.h" "*.cpp" "../../include/sos/*.h" "../../include/sos/*.cpp" "../../include/barriers/*.h" "../../include/barriers/*.cpp"
        "../../external/ChebTools/*"
        "../../external/Padua/*"
        "../../include/sos/utils.cpp"
        )

if (CMAKE_BUILD_TYPE STREQUAL "Release_double_avx_openmp")
    set(USE_OpenMP TRUE)
else ()
    set(USE_OpenMP FALSE)
endif ()

message("USE_OPENMP is ${USE_OpenMP}")

if (APPLE AND USE_OpenMP)
    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        message("libomp library is in ${OpenMP_libomp_LIBRARY}")
    endif ()
endif ()

if (USE_OpenMP)
    find_package(OpenMP REQUIRED)
endif (USE_OpenMP)

find_package(Python2 COMPONENTS Development NumPy)
add_executable(EnvelopeProblem ${TARGETS})

if (NOT BOOST_DIR)
    find_package(Boost 1.67.0 COMPONENTS *boost usr/local/)
    set(BOOST_DIR ${Boost_INCLUDE_DIRS})
    message("Boost found at ${Boost_INCLUDE_DIRS}")
endif ()

if (NOT BOOST_DIR)
    message(FATAL_ERROR "This program requires the boost library, and will not be compiled. Set with flag -BOOST_DIR.")
else ()
    target_include_directories(EnvelopeProblem PRIVATE ${BOOST_DIR}
            ../../include/sos ../../external ${SPDLOG_DIR} ../../include/sos/include ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS}
            ../../../Eigen/eigen
            /usr/local/include
            ${LLVM_INCLUDE_DIRS})
    message("LLVM INCLUDE DIRS are ${LLVM_INCLUDE_DIRS}")
    if (CMAKE_BUILD_TYPE STREQUAL "Release_double_avx_openmp")
        message("Link MKL libraries")
        target_link_options(EnvelopeProblem PRIVATE -framework Accelerate /opt/local/lib/lapack/liblapacke.dylib $ENV{MKLROOT}/lib/libmkl_intel_ilp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
        #Common Error on MacOs is that @rpath/libiomp5.dylib can not be linked/found. Fix is to run "install_name_tool -change @rpath/libiomp5.dylib <MKL_LIB>/lib/libiomp5.dylib ./EXECUTABLE"
    endif ()
    find_package(LLVM REQUIRED CONFIG)
    message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
    message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")
    target_link_directories(EnvelopeProblem PRIVATE ${llvm_libs})
    if (USE_OpenMP)
        target_link_directories(EnvelopeProblem PRIVATE ${OPENMP_LIBRARIES})
    endif ()
    target_link_libraries(EnvelopeProblem Python2::Python Python2::NumPy)
endif ()
