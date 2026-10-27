// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MISC_INTERRUPT_HPP
#define MISC_INTERRUPT_HPP

// Hook for user interrupts (Ctrl-C) in long-running loops. The random walks
// call it once per step. By default it does nothing; a binding defines it
// before including any volesti header, and interrupts by throwing:
//
//   Octave: #define VOLESTI_CHECK_INTERRUPT() octave_quit()
//   R:      #define VOLESTI_CHECK_INTERRUPT() Rcpp::checkUserInterrupt()
//   Python: #define VOLESTI_CHECK_INTERRUPT() \
//               if (PyErr_CheckSignals() != 0) throw pybind11::error_already_set()
//
// It runs once per walk step, so it must be cheap; throttle it in the
// definition if the host's check is not. It is not called from the
// multithreaded walks, where a throw would cross a parallel region.

#ifndef VOLESTI_CHECK_INTERRUPT
#define VOLESTI_CHECK_INTERRUPT() ((void)0)
#endif

#endif // MISC_INTERRUPT_HPP
