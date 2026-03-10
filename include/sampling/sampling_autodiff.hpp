#ifndef SAMPLING_AUTODIFF_H
#define SAMPLING_AUTODIFF_H

// Wrapper templates for autodiff-enabled sampling
// Similar to existing sampling.hpp but optimized for autodiff

#include "sampling.hpp"
#include "ode_solvers/oracle_autodiff_functors.hpp"