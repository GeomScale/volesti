#pragma once
#ifndef BIT_OPERATOR_H
#define BIT_OPERATOR_H
/**
 * @file bit_operator.cpp
 *
 * @brief bit operator utilities
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
 * Algorithms are copied from http://aggregate.org/MAGIC/,
 * collection of public domain algorithms,
 * and modified by authors.
 * You may use original Public Domain Software for your code.
 */
#include "config.h"
#include <inttypes.h>
/**
 * Counts number of 1s.
 *
 * SIMD within a Register algorithm
 * citing from a website http://aggregate.org/MAGIC/
 * @param[in] x bit pattern
 * @return number of 1s in \b x.
 */
static inline int ones(uint32_t x) {
#if HAVE___BUILTIN_POPCOUNTL
    return __builtin_popcountl(x);
#else
    x -= (x >> 1) & UINT32_C(0x55555555);
    x = ((x >> 2) & UINT32_C(0x33333333)) + (x & UINT32_C(0x33333333));
    x = ((x >> 4) + x) & UINT32_C(0x0f0f0f0f);
    x += (x >> 8);
    x += (x >> 16);
    return static_cast<int>(x & UINT32_C(0x3f));
#endif
}

/**
 * Counts number of 1s.
 *
 * SIMD within a Register algorithm
 * citing from a website http://aggregate.org/MAGIC/
 * @param[in] x bit pattern
 * @return number of 1s in \b x.
 */
static inline int ones(uint64_t x) {
#if HAVE___BUILTIN_POPCOUNTLL
    return __builtin_popcountll(x);
#else
    x -= (x >> 1) & UINT64_C(0x5555555555555555);
    x = ((x >> 2) & UINT64_C(0x3333333333333333))
        + (x & UINT64_C(0x3333333333333333));
    x = ((x >> 4) + x) & UINT64_C(0x0f0f0f0f0f0f0f0f);
    x += (x >> 8);
    x += (x >> 16);
    x += (x >> 32);
    return static_cast<int>(x & UINT64_C(0x7f));
#endif
}

#if 0
/**
 * count leading zero from MSB
 * SIMD within a Register algorithm
 * citing from a website http://aggregate.org/MAGIC/
 * @param[in] x bit pattern
 * @return number of 0s from MSB in \b x.
 */
static inline int leadingZeroCount(uint32_t x) {
#if HAVE___BUILTIN_CLZ
    return __builtin_clz(x);
#else
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return 32 - ones(x);
#endif
}

/**
 * count leading zero from MSB
 * SIMD within a Register algorithm
 * citing from a website http://aggregate.org/MAGIC/
 * @param[in] x bit pattern
 * @return number of 0s from MSB in \b x.
 */
static inline int leadingZeroCount(uint64_t x) {
#if HAVE___BUILTIN_CLZL
    return __builtin_clzl(x);
#else
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);
    return 64 - ones(x);
#endif
}
#endif

/**
 * position of least significant one bit
 * SIMD within a Register algorithm
 * citing from a website http://aggregate.org/MAGIC/
 * 00000 -> 0
 * ***10 -> 1
 * **100 -> 2
 * @param[in] x bit pattern
 * @return position of leastSignificantOneBit
 */
static inline int leastSignificantOneBit(uint32_t x) {
    return static_cast<int32_t>(x) & -static_cast<int32_t>(x);
}

/**
 * position of least significant one bit
 * SIMD within a Register algorithm
 * citing from a website http://aggregate.org/MAGIC/
 * 00000 -> 0
 * ***10 -> 1
 * **100 -> 2
 * @param[in] x bit pattern
 * @return position of leastSignificantOneBit
 */
static inline int leastSignificantOneBit(uint64_t x) {
    return static_cast<int64_t>(x) & -static_cast<int64_t>(x);
}

/**
 * position of zero bit from LSB
 * 0000 -> 1
 * **01 -> 2
 * *011 -> 3
 * 1111 -> 5 (should be 0?)
 * @param[in] x bit pattern
 * @return number of 1s from MSB in \b x.
 */
template<typename T>
int firstZeroBit(T x)
{
    return leastSignificantOneBit(~x);
}

template<typename T>
T innerProduct(T a, T b)
{
    return static_cast<T>(ones(a & b) & 1);
}

static inline uint32_t reverseBit(uint32_t x)
{
    x = ((x & UINT32_C(0xaaaaaaaa)) >> 1)
        | ((x & UINT32_C(0x55555555)) << 1);
    x = ((x & UINT32_C(0xcccccccc)) >> 2)
        | ((x & UINT32_C(0x33333333)) << 2);
    x = ((x & UINT32_C(0xf0f0f0f0)) >> 4)
        | ((x & UINT32_C(0x0f0f0f0f)) << 4);
    x = ((x & UINT32_C(0xff00ff00)) >> 8)
        | ((x & UINT32_C(0x00ff00ff)) << 8);
    return (x >> 16) | (x << 16);
}

static inline uint64_t reverseBit(uint64_t x)
{
    x = ((x & UINT64_C(0xaaaaaaaaaaaaaaaa)) >> 1) |
        ((x & UINT64_C(0x5555555555555555)) << 1);
    x = ((x & UINT64_C(0xcccccccccccccccc)) >> 2) |
        ((x & UINT64_C(0x3333333333333333)) << 2);
    x = ((x & UINT64_C(0xf0f0f0f0f0f0f0f0)) >> 4) |
        ((x & UINT64_C(0x0f0f0f0f0f0f0f0f)) << 4);
    x = ((x & UINT64_C(0xff00ff00ff00ff00)) >> 8) |
        ((x & UINT64_C(0x00ff00ff00ff00ff)) << 8);
    x = ((x & UINT64_C(0xffff0000ffff0000)) >> 16) |
        ((x & UINT64_C(0x0000ffff0000ffff)) << 16);
    return (x >> 32) | (x << 32);
}

static inline int getBit(uint64_t x, uint32_t pos)
{
    return (x >> pos) & 1;
}

#endif // BIT_OPERATOR_H
