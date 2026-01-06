// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


// Umbrella header for MCMC and sampling diagnostics.
// Including this header provides access to:
//
//  - Convergence diagnostics (PSRF variants)
//  - Effective sample size estimation
//  - Chain thinning utilities
//  - Diagnostic printing helpers
//
// Intended for:
//  - Post-sampling analysis
//  - MCMC convergence validation
//  - Benchmarking and debugging
//
// Note:
//  - Some diagnostics assume approximately stationary chains
//  - Diagnostics may be computationally expensive for large sample sizes


#ifndef DIAGNOSTICS_DIAGNOSTICS_HPP
#define DIAGNOSTICS_DIAGNOSTICS_HPP

// Formatting utilities used by diagnostic output
#include "misc/print_table.hpp"
// Convergence diagnostics
#include "diagnostics/multivariate_psrf.hpp"
#include "diagnostics/univariate_psrf.hpp"
#include "diagnostics/interval_psrf.hpp"
// Sample quality diagnostics
#include "diagnostics/geweke.hpp"
#include "diagnostics/raftery.hpp"
#include "diagnostics/effective_sample_size.hpp"
// Utilities
#include "diagnostics/thin_samples.hpp"
#include "diagnostics/print_diagnostics.hpp"

#endif
