set(SUITESPARSE_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetSuiteSparse)
    find_path(SUITESPARSE_DIR NAMES include/SuiteSparse_config.h PATHS ${SUITESPARSE_CMAKE_DIR}/../_deps/suitesparse-build/install)

    if (NOT SUITESPARSE_DIR)
        include(FetchContent)
        set(FETCHCONTENT_BASE_DIR "${SUITESPARSE_CMAKE_DIR}/../_deps")
        
        FetchContent_Declare(
            suitesparse
            GIT_REPOSITORY https://github.com/DrTimothyAldenDavis/SuiteSparse.git
            GIT_TAG d558c83006d63d1dc62004f30042b3ca484f3f94
        )

        FetchContent_GetProperties(suitesparse)

        if(NOT suitesparse_POPULATED)
            message(STATUS "SuiteSparse library not found locally, downloading it.")
            FetchContent_Populate(suitesparse)
            
            set(SUITESPARSE_SOURCE_DIR ${suitesparse_SOURCE_DIR})
            set(SUITESPARSE_BUILD_DIR ${suitesparse_BINARY_DIR})
            set(SUITESPARSE_INSTALL_DIR ${SUITESPARSE_BUILD_DIR}/install)

            # Create build directory
            file(MAKE_DIRECTORY ${SUITESPARSE_BUILD_DIR})

            # Configure SuiteSparse with CMake
            # Build only the necessary components: SuiteSparse_config, AMD, COLAMD, CAMD, CCOLAMD, CHOLMOD, and SPQR
            message(STATUS "Configuring SuiteSparse...")
            execute_process(
                COMMAND ${CMAKE_COMMAND} 
                    -S ${SUITESPARSE_SOURCE_DIR}
                    -B ${SUITESPARSE_BUILD_DIR}
                    -DCMAKE_INSTALL_PREFIX=${SUITESPARSE_INSTALL_DIR}
                    -DCMAKE_BUILD_TYPE=Release
                    -DBUILD_SHARED_LIBS=OFF
                    -DBUILD_STATIC_LIBS=ON
                    -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;colamd;camd;ccolamd;cholmod;spqr"
                    -DCHOLMOD_CAMD=ON
                    -DCHOLMOD_SUPERNODAL=ON
                WORKING_DIRECTORY ${SUITESPARSE_BUILD_DIR}
                RESULT_VARIABLE SUITESPARSE_CONFIGURE_RESULT
                OUTPUT_VARIABLE SUITESPARSE_CONFIGURE_OUTPUT
                ERROR_VARIABLE SUITESPARSE_CONFIGURE_ERROR
            )

            if(NOT SUITESPARSE_CONFIGURE_RESULT EQUAL 0)
                message(FATAL_ERROR "SuiteSparse configuration failed: ${SUITESPARSE_CONFIGURE_ERROR}")
            endif()

            # Build SuiteSparse
            message(STATUS "Building SuiteSparse...")
            execute_process(
                COMMAND ${CMAKE_COMMAND} --build ${SUITESPARSE_BUILD_DIR} --config Release
                WORKING_DIRECTORY ${SUITESPARSE_BUILD_DIR}
                RESULT_VARIABLE SUITESPARSE_BUILD_RESULT
                OUTPUT_VARIABLE SUITESPARSE_BUILD_OUTPUT
                ERROR_VARIABLE SUITESPARSE_BUILD_ERROR
            )

            if(NOT SUITESPARSE_BUILD_RESULT EQUAL 0)
                message(FATAL_ERROR "SuiteSparse build failed: ${SUITESPARSE_BUILD_ERROR}")
            endif()

            # Install SuiteSparse to the install directory
            message(STATUS "Installing SuiteSparse...")
            execute_process(
                COMMAND ${CMAKE_COMMAND} --install ${SUITESPARSE_BUILD_DIR} --config Release
                WORKING_DIRECTORY ${SUITESPARSE_BUILD_DIR}
                RESULT_VARIABLE SUITESPARSE_INSTALL_RESULT
                OUTPUT_VARIABLE SUITESPARSE_INSTALL_OUTPUT
                ERROR_VARIABLE SUITESPARSE_INSTALL_ERROR
            )

            if(NOT SUITESPARSE_INSTALL_RESULT EQUAL 0)
                message(FATAL_ERROR "SuiteSparse installation failed: ${SUITESPARSE_INSTALL_ERROR}")
            endif()

            set(SUITESPARSE_DIR ${SUITESPARSE_INSTALL_DIR})
            message(STATUS "SuiteSparse built and installed at: ${SUITESPARSE_DIR}")

        endif()
    else()
        message(STATUS "SuiteSparse library found: ${SUITESPARSE_DIR}")
    endif()

    # Set include directories
    include_directories(BEFORE ${SUITESPARSE_DIR}/include)
    include_directories(BEFORE ${SUITESPARSE_DIR}/include/suitesparse)

    # Set library directory
    set(SUITESPARSE_LIB_DIR ${SUITESPARSE_DIR}/lib PARENT_SCOPE)
    
    # Find the libraries (check for both static and dynamic versions)
    find_library(SPQR_LIB NAMES libspqr.a libspqr.dylib spqr PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    find_library(CHOLMOD_LIB NAMES libcholmod.a libcholmod.dylib cholmod PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    find_library(AMD_LIB NAMES libamd.a libamd.dylib amd PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    find_library(COLAMD_LIB NAMES libcolamd.a libcolamd.dylib colamd PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    find_library(CAMD_LIB NAMES libcamd.a libcamd.dylib camd PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    find_library(CCOLAMD_LIB NAMES libccolamd.a libccolamd.dylib ccolamd PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    find_library(SUITESPARSECONFIG_LIB NAMES libsuitesparseconfig.a libsuitesparseconfig.dylib suitesparseconfig PATHS ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
    
    # Verify all libraries were found
    if(NOT SPQR_LIB OR NOT CHOLMOD_LIB OR NOT AMD_LIB OR NOT COLAMD_LIB OR NOT CAMD_LIB OR NOT CCOLAMD_LIB OR NOT SUITESPARSECONFIG_LIB)
        message(FATAL_ERROR "Failed to find one or more SuiteSparse libraries in ${SUITESPARSE_DIR}/lib")
    else()
        message(STATUS "Found SPQR: ${SPQR_LIB}")
        message(STATUS "Found CHOLMOD: ${CHOLMOD_LIB}")
        message(STATUS "Found AMD: ${AMD_LIB}")
        message(STATUS "Found COLAMD: ${COLAMD_LIB}")
        message(STATUS "Found CAMD: ${CAMD_LIB}")
        message(STATUS "Found CCOLAMD: ${CCOLAMD_LIB}")
        message(STATUS "Found SuiteSparseConfig: ${SUITESPARSECONFIG_LIB}")
    endif()

    # Create an interface library for SuiteSparse
    if(NOT TARGET SuiteSparse::SPQR)
        add_library(SuiteSparse::SPQR INTERFACE IMPORTED GLOBAL)
        target_link_libraries(SuiteSparse::SPQR INTERFACE 
            ${SPQR_LIB}
            ${CHOLMOD_LIB}
            ${AMD_LIB}
            ${COLAMD_LIB}
            ${CAMD_LIB}
            ${CCOLAMD_LIB}
            ${SUITESPARSECONFIG_LIB}
        )
        target_include_directories(SuiteSparse::SPQR INTERFACE 
            ${SUITESPARSE_DIR}/include
            ${SUITESPARSE_DIR}/include/suitesparse
        )
    endif()

    # Export the library list for direct use
    set(SUITESPARSE_LIBRARIES 
        ${SPQR_LIB}
        ${CHOLMOD_LIB}
        ${AMD_LIB}
        ${COLAMD_LIB}
        ${CAMD_LIB}
        ${CCOLAMD_LIB}
        ${SUITESPARSECONFIG_LIB}
        PARENT_SCOPE
    )

    message(STATUS "SuiteSparse libraries configured successfully")

endfunction()

