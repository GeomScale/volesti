struct GaussianHMCCorrelation
{
    GaussianHMCCorrelation(double walk_len, unsigned int _rho)
            :   param(walk_len, true, _rho, true)
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
    typename NT,
    typename RandomNumberGenerator
>
struct Walk
{
    unsigned int _rho;
    NT diam;
    Point _p;
    Point _v;
    NT _lambda_prev;

    Walk(CorreSpectra &P, Point const& p, RandomNumberGenerator &rng,
         parameters const& params)
    {
        diam = P.diam;
        _rho = 100 * P.dimension(); // upper bound for the number of reflections (experimental)
        initialize(P, p, rng);
    }

    inline void apply(CorreSpectra &P,
                      Point& p,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * diam;
            _v = GetDirection<Point>::apply(n, rng, false);
            Point p0 = _p;
            int it = 0;
            while (it < _rho)
            {
                auto pbpair = P.trigonometric_positive_intersect(_p, _v, _omega, _facet_prev);
                if (T <= pbpair.first) {
                    update_position(_p, _v, T, _omega);
                    break;
                }
                _lambda_prev = pbpair.first;
                T -= _lambda_prev;
                update_position(_p, _v, _lambda_prev, _omega);
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == _rho){
                _p = p0;
            }
        }
        p = _p;
    }

private :

    inline void initialize(CorreSpectra & P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        _facet_prev = -1;
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);

        NT T = rng.sample_urdist() * _Len;
        int it = 0;

        while (it <= _rho)
        {
            auto pbpair
                    = P.trigonometric_positive_intersect(_p, _v);
            if (T <= pbpair.first) {
                update_position(_p, _v, T, _omega);
                break;
            } else if (it == _rho) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                update_position(_p, _v, _lambda_prev, _omega);
                break;
            }
            _lambda_prev = pbpair.first;
            update_position(_p, _v, _lambda_prev, _omega);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
    }

    inline void update_position(Point &p, Point &v, NT const& T, NT const& omega)
    {
        NT C, Phi;
        for (size_t i = 0; i < p.dimension(); i++)
        {
            C = std::sqrt(p[i] * p[i] + (v[i] * v[i]) / (omega * omega));
            Phi = std::atan((-v[i]) / (p[i] * omega));
            if (v[i] < 0.0 && Phi < 0.0) {
                Phi += M_PI;
            } else if (v[i] > 0.0 && Phi > 0.0) {
                Phi -= M_PI;
            }
            p.set_coord(i, C * std::cos(omega * T + Phi));
            v.set_coord(i, -C * omega * std::sin(omega * T + Phi));
        }
        
    }
};