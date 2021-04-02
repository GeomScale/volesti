# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis

# Contributed and/or modified by Suraj Choubey

# Licensed under GNU LGPL.3, see LICENCE file


add_executable (vpolytopevolume vpolytopevolume.cpp)
TARGET_LINK_LIBRARIES(vpolytopevolume ${LP_SOLVE})
