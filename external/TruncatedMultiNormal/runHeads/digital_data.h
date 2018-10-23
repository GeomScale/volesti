#pragma once
#ifndef DIGITAL_DATA_H
#define DIGITAL_DATA_H
/**
 * @file digital.h
 *
 * @brief digital net binary file format
 *
 * @author Shinsuke Mori (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Mutsuo Saito (Hiroshima University)
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 */
#include <inttypes.h>

namespace MCQMCIntegration {

    struct dndata {
        uint32_t n;
        uint32_t s;
        uint32_t m;
        uint32_t tvalue;
        double wafom;
        uint64_t data[180]; // max_s * max_m
    };

    extern const dndata nxlw[];
    extern const dndata solw[];
}
#endif // DIGITAL_DATA_H
