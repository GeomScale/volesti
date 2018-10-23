#pragma once
#ifndef POWTWO_H
#define POWTWO_H
/**
 * @file powtwo.h
 *
 * @brief calculate power of two. (2^n)
 *
 * @author Shinsuke Mori (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Mutsuo Saito
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 *
 */

#include <iostream>
#include <iomanip>
#include <inttypes.h>
#include <stdexcept>

static inline uint64_t powtwo(int index) {
    using namespace std;
    if (index < 0 || index > 63) {
        cerr << "index out of range in powtwo. index = " << dec << index
             << "\nindex r should be in the range 0 <= r <= 63" << endl;
        throw invalid_argument("index out of range in powtwo");
    }
    return UINT64_C(1) << index;
}

#endif // POWTWO_H
