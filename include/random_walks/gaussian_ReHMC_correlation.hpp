/// ReHMC random WalkType for sampling correlation matrices

#ifndef RANDOM_WALKS_GAUSSIAN_REHMC_CORRELATION_HPP
#define RANDOM_WALKS_GAUSSIAN_REHMC_CORRELATION_HPP

/// ReHMC walk

struct GaussianReHMCCorrelationWalk
{
    GaussianReHMCCorrelationWalk(double L, unsigned int _rho)
            :   param(L, true, _rho, true)
    {}

    GaussianReHMCCorrelationWalk(double L)
            :   param(L, true, 0, false)
    {}

    GaussianReHMCCorrelationWalk()
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
    typename ConvexBodyType,
    typename RandomNumberGenerator
>
struct Walk
{   
    typedef typename ConvexBodyType::PointType    Point;
    typedef typename Point::FT              NT;

    Walk(ConvexBodyType & P, Point const& p, NT const& a_i, RandomNumberGenerator &rng)
    {   
        _Len = compute_diameter(P);
        _rho = 50 * P.dimension();
        _step = 0.01;
        initialize(P, p, rng);
    }

    Walk(ConvexBodyType & P, Point const& p, NT const& a_i, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                : compute_diameter(P);
        _rho = 50 * P.dimension();
        _step = 0.01;
        initialize(P, p, rng);
    }

    NT compute_diameter(ConvexBodyType const& P){
        std::pair<Point, NT> inner_ball = P.getInnerBall();
        return NT(6) * NT(P.dimension()) * inner_ball.second;
    }

    NT Hamiltonian(Point const &p, Point const &v){
        return (p.dot(p) + v.dot(v))/2;
    }

    void esti_grad(Point const& p, Point & grad_p){
        grad_p = p;
    }

    inline void apply(  ConvexBodyType & P, 
                        Point& p, 
                        NT const& a_i,
                        unsigned int const& walkL, 
                        RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension(), num_leaps = ceil(walkL / _step);
        int i, j, it;
        NT T, h1, h2, u;
        Point p0, grad_p;
        
        for(i = 0; i < walkL; ++i){
            T = _step; // leap total distance T = rng.sample_urdist() * _Len;
            _v = GetDirection<Point>::apply(n, rng, false);
            p0 = _p;

            h1 = Hamiltonian(_p, _v);
            // std::cout << "h1 = " << h1 << std::endl;

            // grad_p = esti_grad(_p, grad_p);

            for(j = 0; j < num_leaps; ++j){            
                _v = _v - (_step/2) * _p;

                // it = 0;
                while(true){
                // while(it < _rho){
                    auto pbpair = P.line_positive_intersect(_p, _v);
                    if (T <= pbpair.first){
                        update_position(_p, _v, T);
                        break;
                    }
                    
                    _lambda_prev = 0.995 * pbpair.first;
                    update_position(_p, _v, _lambda_prev);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, pbpair.second);
                    // it++;
                }
                // if (it == _rho) _p = p0;

                // grad_p = esti_grad(p);
                _v = _v - (_step/2) * _p;
            }
            h2 = Hamiltonian(_p, _v);
            // std::cout << "h2 = " << h2 << std::endl;
            u = ((double)rand()/(double)RAND_MAX);
            if (u < std::exp(h1-h2)) _p = p0;
        }
        p = _p;
    }

    inline void update_position(Point &p, Point &v, NT const& T){
        p += T * v;    
    }

private:

    inline void initialize(ConvexBodyType &P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {   
        unsigned int n = P.dimension();
        
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);
        
        NT T = rng.sample_urdist() * _Len;
        int it = 0;

        while (it <= _rho){
            auto pbpair = P.line_positive_intersect(_p, _v);
            if (T <= pbpair.first){
                _p += T * _v;
                break;
            } else if (it == _rho) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                _p += _lambda_prev * _v;
                break;
            }
            _lambda_prev = 0.995 * pbpair.first;
            _p += _lambda_prev * _v;
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
    }

    NT _Len;
    unsigned int _rho;
    Point _p;
    Point _v;
    NT _lambda_prev;
    NT _step;
    // typename Point::Coeff _lambdas;
    // typename Point::Coeff _Av;
};

};

#endif // RANDOM_WALKS_GAUSSIAN_REHMC_CORRELATION_HPP
