// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef VARS_H
#define VARS_H

//structs with variables and random generators
template <typename NT, typename RNG>
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
          NT diameter,
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
          bool cdhr_walk,
          bool rdhr_walk,
          bool bill_walk
    ) :
            m(m), n(n), walk_steps(walk_steps), n_threads(n_threads), err(err), error(error),
            lw(lw), up(up), L(L), che_rad(che_rad), diameter(diameter), rng(rng),
            urdist(urdist), urdist1(urdist1) , delta(delta) , verbose(verbose), rand_only(rand_only), round(round),
            NN(NN),birk(birk), ball_walk(ball_walk), cdhr_walk(cdhr_walk), rdhr_walk(rdhr_walk), bill_walk(bill_walk){};

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
    NT diameter;
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
    bool cdhr_walk;
    bool rdhr_walk;
    bool bill_walk;
};

template <typename NT, typename RNG>
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
           bool verbose,
           bool rand_only,
           bool round,
           bool NN,
           bool birk,
           bool ball_walk,
           bool cdhr_walk,
           bool rdhr_walk
    ) :
            n(n), walk_steps(walk_steps), N(N), W(W), n_threads(n_threads), error(error),
            che_rad(che_rad), rng(rng), C(C), frac(frac), ratio(ratio), delta(delta),
            verbose(verbose), rand_only(rand_only), round(round),
            NN(NN),birk(birk),ball_walk(ball_walk),cdhr_walk(cdhr_walk), rdhr_walk(rdhr_walk){};

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
    bool verbose;
    bool rand_only;
    bool round;
    bool NN;
    bool birk;
    bool ball_walk;
    bool cdhr_walk;
    bool rdhr_walk;
};


template <typename NT>
struct vars_ban {
public:

    vars_ban(NT lb,
             NT ub,
             NT p,
             NT rmax,
             NT alpha,
             int win_len,
             int N,
             int nu,
             bool window2
    ) :
            lb(lb), ub(ub), p(p), rmax(rmax), alpha(alpha),
            win_len(win_len), N(N), nu(nu), window2(window2) {};


    NT lb;
    NT ub;
    NT p;
    NT rmax;
    NT alpha;
    int win_len;
    int N;
    int nu;
    bool window2;
};



#endif
