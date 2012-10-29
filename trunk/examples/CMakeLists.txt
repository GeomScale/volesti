project( HeaDaCHe_examples )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0005 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories (BEFORE ../external/triangulation/include)
  include_directories (BEFORE ../external/spatial_sorting/include)
  include_directories (BEFORE ../external/kernel_d/include)
  include_directories (BEFORE ../external/eigen)
  include_directories (BEFORE ../include)
  #create_single_source_cgal_program("triangulation.cpp")
  #create_single_source_cgal_program("compare_cgal.cpp")
  
  #include_directories (BEFORE ../patches/include)
  #create_single_source_cgal_program("headache.cpp")
  #create_single_source_cgal_program("eigen_adjoint.cpp")
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11 ")
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")
  #link_libraries("proc")
  create_single_source_cgal_program("opt.cpp")

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
