# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

# Licensed under GNU LGPL.3, see LICENCE file

find_package(OpenMP REQUIRED)

add_executable (skinny_cube_10_dim skinny_cube_10_dim.cpp)
add_executable (random_hpoly_50_dim random_hpoly_50_dim.cpp)
add_executable (random_hpoly_50_dim_parallel random_hpoly_50_dim_parallel.cpp)

TARGET_LINK_LIBRARIES(skinny_cube_10_dim ${LP_SOLVE})
TARGET_LINK_LIBRARIES(random_hpoly_50_dim ${LP_SOLVE})
TARGET_LINK_LIBRARIES(random_hpoly_50_dim_parallel ${LP_SOLVE} OpenMP::OpenMP_CXX)
