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
        _astep = a * _step;
        num_leaps = ceil(_Len / _step);
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
        unsigned int n = P.dimension();
        NT T, h1, h2, u, diff;
        Point p0 = p;

        for(int i = 0; i < walk_length; ++i){
            T = _step;
            _v = GetDirection<Point>::apply(n, rng, false);
            _p = p0;
            h1 = Hamiltonian(_p, _v, a);
            for(int j = 0; j < num_leaps; ++j){   
                _v -= _astep * _p;

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
                _v -= _astep * _p;
            }
            h2 = Hamiltonian(_p, _v, a);
            diff = h1 - h2 < 0 ? h1 - h2 : 0;
            u = log(rng.sample_urdist());
            if(u < diff){
                p0 = _p;
            }
        }
        p = p0;
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }


private :

    inline NT Hamiltonian(Point &p, Point &v, NT a){
        return a * p.squared_length() + 0.5 * v.squared_length();
    }

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
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);
        NT T = _step;
        
        for(int j = 0; j < num_leaps; ++j){            
            _v -= _astep * _p;

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
            _v -= _astep * _p;
        }
        p = _p;
    }

    NT _Len;
    NT _step;
    NT _astep;
    Point _p;
    Point _v;
    NT _lambda_prev;
    unsigned int _rho;
    unsigned int num_leaps;
};

};


#endif // RANDOM_WALKS_GAUSSIAN_REHMC_WALK_HPP

