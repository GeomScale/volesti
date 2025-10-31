// VolEsti (volume computation and sampling library)

// Copyright (c) 2020 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.


#ifndef GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP
#define GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP

#include <chrono>
#include <boost/random.hpp>

namespace detail {
template <typename RNG, typename NT>
inline NT sample_trunc_expdist(RNG& rng, boost::random::exponential_distribution<NT>& expdist) {
    NT z;
    do { z = expdist(rng); } while (z > NT(1));
    return z;
}
}

/////////////////// Random numbers generator
///
/// \tparam RNGType
/// \tparam NT
/// \tparam Ts

template <typename RNGType, typename NT, int ... Ts>
struct BoostRandomNumberGenerator;

template <typename RNGType, typename NT>
struct BoostRandomNumberGenerator<RNGType, NT>
{
    BoostRandomNumberGenerator(int d)
            :   _rng(std::chrono::system_clock::now().time_since_epoch().count())
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
    {}

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

    NT sample_trunc_expdist() 
    {
        return detail::sample_trunc_expdist(_rng, _expdist);
    }

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
    boost::random::exponential_distribution<NT> _expdist;
};


template <typename RNGType, typename NT, int Seed>
struct BoostRandomNumberGenerator<RNGType, NT, Seed>
{
    BoostRandomNumberGenerator(int d=1)
            :   _rng(Seed)
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
    {}

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

    NT sample_trunc_expdist() 
    {
        return detail::sample_trunc_expdist(_rng, _expdist);
    }

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
    boost::random::exponential_distribution<NT> _expdist;
};

#endif // GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP
