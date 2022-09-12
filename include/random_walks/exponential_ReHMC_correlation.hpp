#ifndef RANDOM_WALKS_EXPONENTIAL_REHMC_CORRELATION_HPP
#define RANDOM_WALKS_EXPONENTIAL_REHMC_CORRELATION_HPP

#include "sampling/sphere.hpp"

// Reflective HMC for sampling correlation matrices from the Exponential distribution
// Pdf function is exp(- <c, x+I> / T)

struct ExponentialReHMCWalk
{
    ExponentialReHMCWalk(double L)
            :   param(L, true, 0, false)
    {}

    ExponentialReHMCWalk(double L, unsigned int rho)
            :   param(L, true, rho, true)
    {}

    ExponentialReHMCWalk()
            :   param(0, false, 0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set, unsigned int _rho, bool _set_rho)
                :   m_L(L), set_L(set), rho(_rho), set_rho(_set_rho)
        {}
        double m_L;
        bool set_L;
        unsigned int rho;
        bool set_rho;
    };

    parameters param;

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    template <typename GenericPolytope>
    Walk(GenericPolytope &P, Point const& p, Point const& c, NT const& T, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>::template compute<NT>(P);
        _Temp = T;
        _c = c;
        _step = NT(0.03);
        num_leaps = ceil(_Len / _step);
        _Ac = (_step / (2 * _Temp)) * _c;
        _rho = 100 * P.dimension();
        initialize(P, p, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline bool apply(GenericPolytope &P,
                      Point& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T, h1, h2, u, diff;
        Point p0 = p;

        for(int i = 0; i < walk_length; ++i){
            T = _step;
            _v = GetDirection<Point>::apply(n, rng, false);
            _p = p0;
            h1 = Hamiltonian(_p, _v);

            for(int j = 0; j < num_leaps; ++j){
                _v -= _Ac;

                while(true){
                    auto pbpair = P.line_positive_intersect(_p, _v);

                    if (T <= pbpair.first){
                        _p += T * _v;
                        break;
                    }

                    _lambda_prev = 0.995 * pbpair.first;
                    _p += _lambda_prev * _v;
                    T -= _lambda_prev;

                    P.compute_reflection(_v, _p, pbpair.second);
                }
                _v -= _Ac;
            }
            h2 = Hamiltonian(_p, _v);
            diff = h1 - h2 < 0 ? h1 - h2 : 0;
            u = log(rng.sample_urdist());
            if(u < diff){
                p0 = _p;
            }
        }
        p = p0;
        return true;
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }


private :

    inline NT Hamiltonian(Point &p, Point &v){
        return p.dot(_c)/_Temp + 0.5 * v.squared_length();
    }

    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope & P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension(), it;
        NT T = rng.sample_urdist() * _Len;

        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);

        for(int j = 0; j < num_leaps; ++j){
            _v -= _Ac;

            it = 0;
            while(it <= _rho){
                std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v);

                if (T <= pbpair.first){
                    _p += T * _v;
                    _lambda_prev = T;
                    break;
                }

                _lambda_prev = 0.995 * pbpair.first;
                _p += _lambda_prev * _v;
                T -= _lambda_prev;

                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            _v -= _Ac;
        }
    }

    NT _Len;
    NT _step;
    Point _p;
    Point _v;
    Point _c;
    Point _Ac;
    NT _Temp;
    NT _lambda_prev;
    unsigned int _rho;
    unsigned int num_leaps;
};

};


#endif // RANDOM_WALKS_EXPONENTIAL_REHMC_CORRELATION_HPP

