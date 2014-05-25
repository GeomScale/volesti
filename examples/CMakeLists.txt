project( GeomRandWalks_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

find_package(CGAL QUIET COMPONENTS Core )

#FIND_PACKAGE( OpenMP REQUIRED)
#if(OPENMP_FOUND)
#message("OPENMP FOUND")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  #include_directories (BEFORE ../external/extreme_points_d/include)
  include_directories (BEFORE ../external/eigen_3_2)
  include_directories (BEFORE ../include)

  
  #include_directories (BEFORE ../patches/include)
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")
  #link_libraries("proc")
  #create_single_source_cgal_program("opt.cpp")
  create_single_source_cgal_program("vol.cpp")
  #create_single_source_cgal_program("min_ellipsoid.cpp")
  #create_single_source_cgal_program("polytope_generators/cube.cpp")
  #create_single_source_cgal_program("polytope_generators/cross.cpp")
  #create_single_source_cgal_program("polytope_generators/skinny_cube.cpp")
  #create_single_source_cgal_program("polytope_generators/simplex.cpp")
  create_single_source_cgal_program("polytope_generators/birk.cpp")
  #create_single_source_cgal_program("NN_cgal.cpp")
  #create_single_source_cgal_program("mink_sum_vol.cpp")

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
