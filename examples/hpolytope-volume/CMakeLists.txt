# VolEsti (volume computation and sampling library)

# Copyright (c) 20012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis

# Contributed and/or modified by Vaibhav Thakkar

# Licensed under GNU LGPL.3, see LICENCE file


add_executable (hpolytopeVolume hpolytopeVolume.cpp)
TARGET_LINK_LIBRARIES(hpolytopeVolume ${LP_SOLVE})