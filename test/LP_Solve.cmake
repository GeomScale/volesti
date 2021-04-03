function(GetLPSolve)
    find_path(LP_SOLVE_DIR NAMES LPsolve_src PATHS ../../external)

    if (NOT LP_SOLVE_DIR)    
        include(FetchContent)
        FetchContent_Declare(
            lpsolve 
            URL https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz/download 
            URL_HASH MD5=a829a8d9c60ff81dc72ff52363703886
        )
        
        FetchContent_GetProperties(lpsolve)
            
        if(NOT lpsolve_POPULATED)
            message(STATUS "lp_solve library not found locally, downloading it.")
            FetchContent_Populate(lpsolve)
        endif()
        
        set(LP_SOLVE_DIR "${lpsolve_SOURCE_DIR}")
        message(STATUS "Using downloaded lp_solve at: ${LP_SOLVE_DIR}")

    else()

        message(STATUS "lp_solve library found: ${LP_SOLVE_DIR}")
    
    endif()

    set(ENV{LP_SOLVE_DIR})
        
endfunction()
