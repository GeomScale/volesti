#ifndef RANDOM_WALKS_GAUSSIAN_REHMC_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_REHMC_WALK_HPP

#include "sampling/sphere.hpp"

// Reflective HMC for sampling correlation matrices from the Gaussian distribution
// Pdf function is exp(- a * |x|^2)

struct GaussianReHMCWalk
{
    GaussianReHMCWalk(double L)
            :   param(L, true, 0, false)
    {}

    GaussianReHMCWalk(double L, unsigned int rho)
            :   param(L, true, rho, true)
    {}

    GaussianReHMCWalk()
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
    Walk(GenericPolytope &P, Point & p, NT const& a, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>::template compute<NT>(P);
        _step = NT(0.05);
        _rho = 100 * P.dimension();
        initialize(P, p, a, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope &P,
                      Point& p,   // a point to start
                      NT a,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension(), num_leaps = ceil(_Len / _step);
        NT T, h1, h2, u, diff;
        Point p0 = p;

        h1 = Hamiltonian(_p, _v, a);

        for(int i = 0; i < walk_length; ++i){
            T = _step;
            _v = GetDirection<Point>::apply(n, rng, false);
            _p = p0;
            
            for(int j = 0; j < num_leaps; ++j){   
                _v -= (_step * a) * _p;

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
                _v -= (_step * a) * _p;
            }
            h2 = Hamiltonian(_p, _v, a);
            diff = h1 - h2 < 0 ? h1 - h2 : 0;
            u = log(rng.sample_urdist());
            if(u < diff){
                p0 = _p;
                h1 = h2;
            }
        }
        p = p0;
    }

    inline NT Hamiltonian(Point &p, Point &v, NT a){
        return a * p.squared_length() + 0.5 * v.squared_length();
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }


private :

    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope & P,
                           Point & p,
                           NT a,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension(), it;
        NT T = rng.sample_urdist() * _Len;
        unsigned int num_leaps = ceil(T / _step);

        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);

        for(int j = 0; j < num_leaps; ++j){            
            _v -= (_step * a) * _p;

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
            _v -= (_step * a) * _p;
        }
        p = _p;
    }

    NT _Len;
    NT _step;
    Point _p;
    Point _v;
    NT _lambda_prev;
    unsigned int _rho;
};

};


#endif // RANDOM_WALKS_GAUSSIAN_REHMC_WALK_HPP

