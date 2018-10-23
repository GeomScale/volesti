#ifndef KAHAN_HPP
#define KAHAN_HPP
/**
 * @file kahan.hpp
 *
 * @brief Kahan summation algorithm
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
 * @see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */

/*
 * Kahan summation algorithm
 */
class Kahan {
public:
    Kahan() {
        sum = 0.0;
        c = 0.0;
    }
    void clear() {
        sum = 0.0;
        c = 0.0;
    }
    void add(const double x) {
        const double y = x - c;
        const double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    double get() {
        return sum;
    }
private:
    double c;
    double sum;
};

#endif // KAHAN_HPP
