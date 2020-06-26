// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

#include "euler.hpp"
#include "runge_kutta.hpp"
#include "leapfrog.hpp"
#include "bulirsch_stoer.hpp"

#ifndef DISABLE_NLP_ORACLES
#include "collocation.hpp"
#endif

#ifndef ODE_SOLVERS_HPP
#define ODE_SOLVERS_HPP

#endif
