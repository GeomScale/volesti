#ifndef RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <stdexcept>
#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/convex_body.h"
#include "random_walks/compute_diameter.hpp"
#include "random_walks/accelerated_billiard_walk_utils.hpp"

struct BilliardShakeAndBakeWalk
{
    BilliardShakeAndBakeWalk(double L) : param(L, true) {}
    BilliardShakeAndBakeWalk() : param(0, false) {}

    struct parameters
    {
        parameters(double L, bool set) : m_L(L), set_L(set) {}
        double m_L;  bool set_L;
    };

    parameters param;

    template <typename Polytope, typename RandomNumberGenerator>
    struct Walk
    {
        using Point    = typename Polytope::PointType;
        using NT       = typename Point::FT;
        using VT       = typename Polytope::VT;
        using MT       = typename Polytope::MT;
        using MTdense  = Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>;

        struct update_parameters { int facet_prev{-1}; } params_;
        static constexpr NT kDefaultEps = NT(1e-10);

        template <typename GenericPolytope>
        Walk(GenericPolytope &P, const Point &boundary_pt, int facet_idx,
             RandomNumberGenerator &rng, NT eps = kDefaultEps)
        : P_(P), epsilon_(eps)
        {
            if(!P_.is_normalized()) P_.normalize();
            _Len = compute_diameter<GenericPolytope>::template compute<NT>(P_);
            initialize(boundary_pt, facet_idx);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope &P, const Point &boundary_pt, int facet_idx,
             RandomNumberGenerator &rng, NT eps, parameters const& par)
        : P_(P), epsilon_(eps)
        {
            if(!P_.is_normalized()) P_.normalize();
            _Len = par.set_L ? par.m_L : compute_diameter<GenericPolytope>::template compute<NT>(P_);
            initialize(boundary_pt, facet_idx);
        }

        NT get_epsilon() const noexcept { return epsilon_; }


        void apply(unsigned walk_len, RandomNumberGenerator& rng)
        {
            for (unsigned step = 0; step < walk_len; ++step)
            {
                params_.facet_prev = -1;          
                _lambda_prev       = NT(0);       
                lambda_hit_        = NT(0);       


                Point v = get_direction(rng);
                NT T = rng.sample_urdist() * _Len;
                //THIS KEEP STAYING ZERO :(
                Point p0 = p_;
                int reflections = 0;

                Av_ = P_.get_mat() * v.getCoefficients();
                // I know I shouldnt compute this like this :)
                Ar_ = P_.get_mat() * p_.getCoefficients();

                while (T > epsilon_ && reflections < 50*dim_)
                {
                    NT  lambda_min;
                    int facet_new;
                    int attempts = 0;

                    do {
                        std::tie(lambda_min, facet_new) =
                            P_.line_positive_intersect(p_, v, Ar_, Av_, lambda_hit_);
                    } while ((lambda_min <= epsilon_ ||
                            facet_new == params_.facet_prev)   
                            && ++attempts < dim_);

                    if (lambda_min <= epsilon_ || facet_new < 0) break;

                    if (T <= lambda_min) {                     
                        p_  += T * v;
                        Ar_ += T * Av_;                        
                        T    = NT(0);
                    } else {                                   
                        p_  += lambda_min * v;
                        Ar_ += lambda_min * Av_;
                        T   -= lambda_min;

                        P_.compute_reflection(v, p_, facet_new); 
                        Av_ = P_.get_mat() * v.getCoefficients();

                        params_.facet_prev = facet_new;        
                        lambda_hit_        = NT(0);            
                        _lambda_prev       = lambda_min;       
                        reflections++;
                    }
                }

                if (reflections >= 50*dim_)  p_ = p0;         
            }

        }


        const Point& getCurrentPoint() const noexcept { return p_; }
        void update_delta(NT L){ _Len=L; }

    private:

        Point get_direction(RandomNumberGenerator& rng)
        {
            VT z = GetDirection<Point>::apply(dim_, rng).getCoefficients();
            MT I_cc = - A_row_k_ * A_row_k_.transpose();
            I_cc.diagonal() += VT::Ones(dim_);
            NT U = rng.sample_urdist();               
            NT r = std::pow(U, NT(1)/NT(dim_-1)); 
            NT cz = A_row_k_.dot(z);
            VT z_tilde  = I_cc*z;
            z_tilde *= r;
            z_tilde /= std::sqrt(NT(1) - cz*cz);
            
            VT v = z_tilde - std::sqrt(NT(1) - r*r) * A_row_k_;
            return Point(v);
        }


        void initialize(const Point& boundary_pt,
                        int   facet_idx,
                        RandomNumberGenerator& rng)
        {
            dim_ = P_.dimension();
            m_ = P_.num_of_hyperplanes();
            VT b=P_.get_vec();

            NT kFacetEps = epsilon_;

            // Checking if facet index belongs to the boundary point 
            p_ = boundary_pt;
            VT ai = P_.get_facet_normal_vec(facet_idx);
            NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(facet_idx));
            if (dist > kFacetEps)
            {
                facet_idx_ = -1;
                for (int i = 0; i < m_; ++i) {
                    VT ai = P_.get_facet_normal_vec(i);
                    NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(i));
                    if (dist < kFacetEps) {
                        facet_idx_ = i;
                        break;
                    }
                }
                if (facet_idx_ < 0)
                {
                    throw std::runtime_error("Boundary point not on any facet!");
                }
            }
            facet_idx_ = facet_idx;

            //Normal of active facet
            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);

            //Calculating first Ar and initializing Av 
            Ar_.setZero(m_);
            Av_.setZero(m_);
            lambda_hit_ = NT(0);
            
            Ar_.noalias() = P_.get_mat() * p_.getCoefficients();
            lambda_hit_ = NT(0);

            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);

            params_.facet_prev  = facet_idx_;

        }

        Polytope &P_;
        NT epsilon_{kDefaultEps};
        int dim_{0}, m_{0};
        NT    _lambda_prev = NT(0);  
        Point p_;                
        VT    Ar_, Av_;           
        VT    A_row_k_;           
        NT    lambda_hit_{0};
        NT    _Len{1};
    };
};

#endif /* RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP */