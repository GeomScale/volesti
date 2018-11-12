// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file
#ifndef VARS_H
#define VARS_H

//An implementation of Welford's algorithm for mean and variance.
template <typename NT>
std::pair<NT, NT> getMeanVariance(std::vector<NT>& vec) {
    NT mean = 0, M2 = 0, variance = 0, delta;
    typedef typename std::vector<NT>::iterator viterator;

    unsigned int i=0;
    viterator vecit = vec.begin();
    for( ; vecit!=vec.end(); vecit++, i++){
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }

    return std::pair<NT, NT> (mean, variance);
}


//structs with variables and random generators
template <typename NT, class RNG>
struct vars{
public:
    typedef RNG RNGType;
    vars( unsigned int m,
          unsigned int n,
          unsigned int walk_steps,
          unsigned int n_threads,
          const NT err,
          NT error,
          const int lw,
          NT up,
          const int L,
          NT che_rad,
          RNG &rng,
          boost::random::uniform_real_distribution<>(urdist),
          boost::random::uniform_real_distribution<> urdist1,
          NT delta,
          bool verbose,
          bool rand_only,
          bool round,
          bool NN,
          bool birk,
          bool ball_walk,
          bool coordinate
    ) :
            m(m), n(n), walk_steps(walk_steps), n_threads(n_threads), err(err), error(error),
            lw(lw), up(up), L(L), che_rad(che_rad), rng(rng),
            urdist(urdist), urdist1(urdist1) , delta(delta) , verbose(verbose), rand_only(rand_only), round(round),
            NN(NN),birk(birk), ball_walk(ball_walk), coordinate(coordinate){};

    unsigned int m;
    unsigned int n;
    unsigned int walk_steps;
    unsigned int n_threads;
    const NT err;
    NT error;
    const int lw;
    NT up;
    const int L;
    NT che_rad;
    RNG &rng;
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1;
    NT delta;
    bool verbose;
    bool rand_only;
    bool round;
    bool NN;
    bool birk;
    bool ball_walk;
    bool coordinate;
};

template <typename NT, class RNG>
struct vars_g{
public:
    typedef RNG RNGType;
    vars_g(unsigned int n,
           unsigned int walk_steps,
           unsigned int N,
           unsigned int W,
           unsigned int n_threads,
           NT error,
           NT che_rad,
           RNG &rng,
           NT C,
           NT frac,
           NT ratio,
           NT delta,
           bool deltaset,
           bool verbose,
           bool rand_only,
           bool round,
           bool NN,
           bool birk,
           bool ball_walk,
           bool coordinate
    ) :
            n(n), walk_steps(walk_steps), N(N), W(W), n_threads(n_threads), error(error),
            che_rad(che_rad), rng(rng), C(C), frac(frac), ratio(ratio), delta(delta),
            deltaset(deltaset), verbose(verbose), rand_only(rand_only), round(round),
            NN(NN),birk(birk),ball_walk(ball_walk),coordinate(coordinate){};

    unsigned int n;
    unsigned int walk_steps;
    unsigned int N;
    unsigned int W;
    unsigned int n_threads;
    NT error;
    NT che_rad;
    RNG &rng;
    NT C;
    NT frac;
    NT ratio;
    NT delta;
    bool deltaset;
    bool verbose;
    bool rand_only;
    bool round;
    bool NN;
    bool birk;
    bool ball_walk;
    bool coordinate;
};

#endif
