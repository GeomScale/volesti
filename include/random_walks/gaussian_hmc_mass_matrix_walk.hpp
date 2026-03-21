// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// Exact HMC for sampling from Gaussian distribution with a general mass matrix.
// Based on:
//   [10] Chevallier, Pion, Cazals - "Hamiltonian Monte Carlo with boundary reflections,
//        and application to polytope volume calculations"
//   [11] Pakman, Paninski - "Exact Hamiltonian Monte Carlo for Truncated Multivariate Gaussians"

#ifndef RANDOM_WALKS_GAUSSIAN_HMC_MASS_MATRIX_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_HMC_MASS_MATRIX_WALK_HPP

#include "sampling/sphere.hpp"
#include <Eigen/Dense>

// Exact HMC for sampling from the Gaussian distribution N(0, Sigma)
// truncated to a convex polytope, with a general mass matrix M.
//
// The Hamiltonian is:
//   H(x, p) = 0.5 * x^T * Sigma^{-1} * x + 0.5 * p^T * M^{-1} * p
//
// When M = I, this reduces to the standard GaussianHamiltonianMonteCarloExactWalk.
// For order polytopes, M can be chosen to exploit the DAG structure.

struct GaussianHMCMassMatrixWalk
{
    GaussianHMCMassMatrixWalk(double L, unsigned int _rho)
            :   param(L, true, _rho, true)
    {}

    GaussianHMCMassMatrixWalk(double L)
            :   param(L, true, 0, false)
    {}

    GaussianHMCMassMatrixWalk()
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
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    template <typename GenericPolytope>
    Walk(GenericPolytope &P, Point const& p, NT const& a_i, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        _rho = 100 * P.dimension(); // upper bound for reflections
        _use_mass_matrix = false;
        initialize(P, p, a_i, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope &P, Point const& p, NT const& a_i, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        _rho = params.set_rho ? params.rho : 100 * P.dimension();
        _use_mass_matrix = false;
        initialize(P, p, a_i, rng);
    }

    // Constructor with mass matrix M.
    // Momentum is drawn from N(0, M), i.e. p = L_chol * z, z ~ N(0,I)
    // where M = L_chol * L_chol^T (Cholesky decomposition).
    template <typename GenericPolytope>
    Walk(GenericPolytope &P, Point const& p, NT const& a_i, RandomNumberGenerator &rng,
         parameters const& params, MT const& mass_matrix)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        _rho = params.set_rho ? params.rho : 100 * P.dimension();
        _use_mass_matrix = true;

        // Cholesky decomposition of the mass matrix
        Eigen::LLT<MT> llt(mass_matrix);
        _L_chol = llt.matrixL();
        _M_inv = mass_matrix.inverse();

        initialize(P, p, a_i, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      Point& p,
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T;

        GenericPolytope P_normalized = P;
        P_normalized.normalize();

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;

            // Generate momentum: if mass matrix is set, p = L*z; else p ~ N(0,I)
            if (_use_mass_matrix) {
                VT z(n);
                for (unsigned int k = 0; k < n; ++k) {
                    z(k) = rng.sample_ndist();
                }
                VT momentum = _L_chol * z;
                _v = Point(momentum);
            } else {
                _v = GetDirection<Point>::apply(n, rng, false);
            }

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
                nudge_in(P_normalized, _p);
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == _rho){
                _p = p0;
            }
        }
        p = _p;
    }


    template
    <
        typename GenericPolytope
    >
    inline void get_starting_point(GenericPolytope const& P,
                                   Point const& center,
                                   Point &q,
                                   unsigned int const& walk_length,
                                   RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT radius = P.InnerBall().second;

        q = GetPointInDsphere<Point>::apply(n, radius, rng);
        q += center;
        initialize(P, q, rng);

        apply(P, q, walk_length, rng);
    }


    template
    <
        typename GenericPolytope
    >
    inline void parameters_burnin(GenericPolytope const& P,
                                  Point const& center,
                                  unsigned int const& num_points,
                                  unsigned int const& walk_length,
                                  RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        Point p(n);
        std::vector<Point> pointset;
        pointset.push_back(center);
        pointset.push_back(_p);
        NT rad = NT(0), max_dist, Lmax = NT(0), radius = P.InnerBall().second;

        for (int i = 0; i < num_points; i++)
        {
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
            p += center;
            initialize(P, p, rng);

            apply(P, p, walk_length, rng);
            max_dist = get_max_distance(pointset, p, rad);
            if (max_dist > Lmax)
            {
                Lmax = max_dist;
            }
            if (2.0*rad > Lmax)
            {
                Lmax = 2.0 * rad;
            }
            pointset.push_back(p);
        }

        if (Lmax > _Len)
        {
            update_delta(Lmax);
        }
        pointset.clear();
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
    inline void initialize(GenericPolytope const& P,
                           Point const& p,
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        _facet_prev = -1;
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);

        NT T = rng.sample_urdist() * _Len;
        int it = 0;

        GenericPolytope P_normalized = P;
        P_normalized.normalize();

        while (it <= _rho)
        {
            auto pbpair
                    = P.trigonometric_positive_intersect(_p, _v, _omega, _facet_prev);
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
            nudge_in(P_normalized, _p);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
    }

    template
    <
        typename GenericPolytope
    >
    inline void nudge_in(GenericPolytope& P, Point& p, NT tol=NT(0))
    {
        MT A_mat = P.get_dense_mat();
        VT b_vec = P.get_vec();
        int m = A_mat.rows();

        VT b_Ax = b_vec - A_mat * p.getCoefficients();
        const NT* b_Ax_data = b_Ax.data();

        NT dist;

        for (int i = 0; i < m; i++) {
            dist = *b_Ax_data;
            if (dist < NT(-tol)){
                NT eps = -1e-7;
                NT eps_1 = -dist;
                VT A_i = A_mat.row(i);
                NT eps_2 = eps_1 + eps;
                Point shift(A_i);
                shift.operator*=(eps_2);
                p.operator+=(shift);
            }
            b_Ax_data++;
        }
    }

    inline void update_position(Point &p, Point &v, NT const& T, NT const& omega)
    {
        NT next_p, next_v;

        NT sinVal = std::sin(omega * T);
        NT cosVal = std::cos(omega * T);

        NT factor1 = sinVal / omega;
        NT factor2 = -omega * sinVal;

        for (size_t i = 0; i < p.dimension(); i++)
        {
            next_p = cosVal * p[i] + v[i] * factor1;
            next_v = factor2 * p[i] + cosVal * v[i];

            p.set_coord(i, next_p);
            v.set_coord(i, next_v);
        }
    }

    inline double get_max_distance(std::vector<Point> &pointset, Point const& q, double &rad)
    {
        double dis = -1.0, temp_dis;
        int jj = 0;
        for (auto vecit = pointset.begin(); vecit!=pointset.end(); vecit++, jj++)
        {
            temp_dis = (q.getCoefficients() - (*vecit).getCoefficients()).norm();
            if (temp_dis > dis) {
                dis = temp_dis;
            }
            if (jj == 0 && temp_dis > rad) {
                rad = temp_dis;
            }
        }
        return dis;
    }

    int _facet_prev;
    unsigned int _rho;
    NT _Len;
    Point _p;
    Point _v;
    NT _omega;
    NT _lambda_prev;
    bool _use_mass_matrix;
    MT _L_chol;   // Cholesky factor of mass matrix
    MT _M_inv;    // Inverse of mass matrix
};

};


#endif // RANDOM_WALKS_GAUSSIAN_HMC_MASS_MATRIX_WALK_HPP
