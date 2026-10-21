// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef LP_ORACLE_OPTIONS_HPP
#define LP_ORACLE_OPTIONS_HPP

#include "Highs.h"
#include <functional>

// Callback used for the user to apply his own settings. When one
// is given the oracles apply it instead of the defaults, so the user is
// responsible for every option he cares about.
//
// Note: oracles need a vertex solution, so be cautious when using solvers
// other than simplex.
using LPOracleOptions = std::function<void(Highs&)>;

// Configures highs for the lp oracles.
inline void lp_oracles_configure_highs(Highs & highs, LPOracleOptions const& opts = nullptr) {
    if (opts) {
        opts(highs);
    } else {
        highs.setOptionValue("output_flag", false);
        highs.setOptionValue("solver", "simplex");
    }
}

// The result of an lp oracle.
//
// value is meaningful only when solved is true, so callers should
// test the result before reading it.
// @tparam T the type of the value the oracle computes
template <typename T>
struct LPOracleResult {
    T value{};
    bool solved = false;

    explicit operator bool() const noexcept {return solved;}
};

#endif