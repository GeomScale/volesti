// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_ACCELERATED_SPEEDUP_BILLIARD_WALK_HPP
#define RANDOM_WALKS_ACCELERATED_SPEEDUP_BILLIARD_WALK_HPP

#include "convex_bodies/hpolytope.h"


// Billiard walk which accelarates each step for uniform distribution

struct AcceleratedSpeedpBilliardWalk
{
    AcceleratedSpeedpBilliardWalk(double L)
            :   param(L, true)
    {}

    AcceleratedSpeedpBilliardWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
    };

    parameters param;


    template
            <
                    typename Polytope,
                    typename RandomNumberGenerator
            >
    struct Walk
    {
        //typedef typename Polytope::PointType Point;
        typedef typename Polytope::MT MT;
        typedef typename Polytope::VT VT;
        typedef typename Polytope::NT NT;

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, VT const& p, RandomNumberGenerator &rng)
        {
            _update_parameters = update_parameters();
            _L = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
            //std::cout<<"computing AA..."<<std::endl;
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            //std::cout<<"AA computed!"<<std::endl;
            initialize(P, p, rng);
            //std::cout<<"Initialization of the random walk done!"<<std::endl;
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, VT const& p, RandomNumberGenerator &rng,
             parameters const& params)
        {
            _update_parameters = update_parameters();
            _L = params.set_L ? params.m_L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            //std::cout<<"computing AA..."<<std::endl;
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            //std::cout<<"AA computed!"<<std::endl;
            initialize(P, p, rng);
            //std::cout<<"Initialization of the random walk done!"<<std::endl;
        }

        template
                <
                        typename GenericPolytope
                >
        inline void apply(GenericPolytope const& P,
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T;//, avg_ref = 0.0;
            const NT dl = 0.995;
            int it;

            for (auto j=0u; j<walk_length; ++j)
            {
                T = -std::log(rng.sample_urdist()) * _L;
                GetDirectionVT<VT>::apply(_v, rng);
                _p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;

                while (it < 250*n)
                {
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;
                }
                //avg_ref += NT(it);
                if (it == 250*n){
                    //std::cout<<"reflection limit reached"<<std::endl;
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        //break;
                    } else {
                        _lambda_prev = rng.sample_urdist() * pbpair.first;
                        _p += (_lambda_prev * _v);
                    }
                    //_p = _p0;
                } 
            }
            //avg_ref *= (1.0 / (NT(walk_length)));
            //std::cout<<"avg_ref = "<<avg_ref<<std::endl;
            //p = _p;
        }

        inline double get_max_distance(std::vector<VT> &pointset, VT const& q, double &rad) 
        {
            double dis = -1.0, temp_dis;
            int jj = 0;
            for (auto vecit = pointset.begin(); vecit!=pointset.end(); vecit++, jj++) 
            {
                temp_dis = (q - (*vecit)).norm();
                if (temp_dis > dis) {
                    dis = temp_dis;
                }
                if (jj == 0) {
                    if (temp_dis > rad) {
                        rad = temp_dis;
                    }
                }
            }
            return dis;
        }

        template
                <
                        typename GenericPolytope
                >
        inline void get_starting_point(GenericPolytope const& P,
                           VT const& center,
                           VT &q,
                           unsigned int const& walk_length,
                           RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T, Lmax = _L, max_dist, rad = 0.0, radius = P.InnerBall().second;
            const NT dl = 0.995;
            int it;
            //std::vector<VT> pointset;
            /// pointset.push_back(center);
            //pointset.push_back(_p);
            //_p0 = _p;

            for (int i = 0; i < 1; i++) 
            {
                VT p = GetPointInDsphereVT<VT>::apply(n, radius, rng);
                p += center;
                initialize(P, p, rng);
                for (auto j=0u; j<walk_length; ++j)
                {
                    //update_delta(Lmax);
                    T = -std::log(rng.sample_urdist()) * _L;
                    GetDirectionVT<VT>::apply(_v, rng);
                    _p0 = _p;

                    it = 0;
                    std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                                          _update_parameters);
                    
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        continue;
                    }

                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;

                    while (it < 250*n)
                    {
                        std::pair<NT, int> pbpair
                                = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                            _AA, _update_parameters);
                        if (T <= pbpair.first) {
                            _p += (T * _v);
                            _lambda_prev = T;
                            break;
                        }
                        _lambda_prev = dl * pbpair.first;
                        _p += (_lambda_prev * _v);
                        T -= _lambda_prev;
                        P.compute_reflection(_v, _p, _update_parameters);
                        it++;
                    }
                    if (it == 250*n){
                        //std::cout<<"reflection limit reached"<<std::endl;
                        std::pair<NT, int> pbpair
                                = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                            _AA, _update_parameters);
                        if (T <= pbpair.first) {
                            _p += (T * _v);
                            _lambda_prev = T;
                            //break;
                        } else {
                            _lambda_prev = rng.sample_urdist() * pbpair.first;
                            _p += (_lambda_prev * _v);
                        }
                        //_p = _p0;
                    } 
                }
                //std::cout<<"[starting_point], _p is_in = "<<P.is_in(_p)<<std::endl;
                //p = _p;
            }
            //std::cout<<"get q"<<std::endl;
            q = _p;
        }

        template
                <
                        typename GenericPolytope
                >
        inline void burnin(GenericPolytope const& P,
                           VT const& center,
                           unsigned int const& num_points,
                           unsigned int const& walk_length,
                           NT &avg_ref,
                           RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T, Lmax = _L, max_dist, rad = 0.0, radius = P.InnerBall().second;
            const NT dl = 0.995;
            int it;
            std::vector<VT> pointset;
            pointset.push_back(center);
            pointset.push_back(_p);
            //_p0 = _p;

            for (int i = 0; i < num_points; i++) 
            {
                VT p = GetPointInDsphereVT<VT>::apply(n, radius, rng);
                p += center;
                initialize(P, p, rng);
                for (auto j=0u; j<walk_length; ++j)
                {
                    //update_delta(Lmax);
                    T = -std::log(rng.sample_urdist()) * _L;
                    GetDirectionVT<VT>::apply(_v, rng);
                    _p0 = _p;

                    it = 0;
                    std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                                          _update_parameters);
                    if (pbpair.first > Lmax) {
                        Lmax = pbpair.first;
                        //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                    }
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        max_dist = get_max_distance(pointset, _p, rad);
                        if (max_dist > Lmax) {
                            Lmax = max_dist;
                            //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                        }
                        if (2.0*rad > Lmax) {
                            Lmax = 2.0*rad;
                            //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                        }
                        pointset.push_back(_p);
                        _lambda_prev = T;
                        continue;
                    }

                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;

                    while (it < 250*n)
                    {
                        std::pair<NT, int> pbpair
                                = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                            _AA, _update_parameters);
                        if (pbpair.first > Lmax) {
                            Lmax = pbpair.first;
                            //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                            //update_delta(Lmax);
                        }
                        if (T <= pbpair.first) {
                            _p += (T * _v);
                            _lambda_prev = T;
                            max_dist = get_max_distance(pointset, _p, rad);
                            if (max_dist > Lmax) {
                                Lmax = max_dist;
                                //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                            }
                            if (2.0*rad > Lmax) {
                                Lmax = 2.0*rad;
                                //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                            }
                            pointset.push_back(_p);
                            break;
                        }
                        _lambda_prev = dl * pbpair.first;
                        _p += (_lambda_prev * _v);
                        T -= _lambda_prev;
                        P.compute_reflection(_v, _p, _update_parameters);
                        it++;
                    }
                    avg_ref += NT(it);
                    if (it == 250*n){
                        //std::cout<<"reflection limit reached"<<std::endl;
                        std::pair<NT, int> pbpair
                                = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                            _AA, _update_parameters);
                        if (pbpair.first > Lmax) {
                            Lmax = pbpair.first;
                            //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                            //update_delta(Lmax);
                        }
                        if (T <= pbpair.first) {
                            _p += (T * _v);
                            _lambda_prev = T;
                            //break;
                        } else {
                            _lambda_prev = rng.sample_urdist() * pbpair.first;
                            _p += (_lambda_prev * _v);
                        }
                        max_dist = get_max_distance(pointset, _p, rad);
                        if (max_dist > Lmax) {
                            Lmax = max_dist;
                            //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                        }
                        if (2.0*rad > Lmax) {
                            Lmax = 2.0*rad;
                            //std::cout<<"L updated, Lmax = "<<Lmax<<std::endl;
                        }
                        pointset.push_back(_p);
                        //_p = _p0;
                    } 
                }
                //starting_points.col(i) = _p;
                //std::cout<<"[burnin], _p is_in = "<<P.is_in(_p)<<std::endl;
                //p = _p;
            }
            avg_ref *= (1.0 / (NT(num_points) * NT(walk_length)));
            if (Lmax > _L) {
                //std::cout<<"we now update L with the velue: "<<Lmax<<std::endl;
                if (n <= 2500) {
                    update_delta(Lmax);
                }
            }
            pointset.clear();
        }

        inline void update_delta(NT L)
        {
            _L = L;
        }

        inline VT get_curr_sample() {
            return _p;
        }

    

        template
                <
                        typename GenericPolytope
                >
        inline void initialize(GenericPolytope const& P,
                               VT const& p,
                               RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _p0.setZero(n);
            _v.setZero(n);
            GetDirectionVT<VT>::apply(_v, rng);

            NT T = -std::log(rng.sample_urdist()) * _L;
            //p0 = _p;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, _update_parameters);

            while (it <= 200*n)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                } else if (it == 200*n) {
                    _lambda_prev = rng.sample_urdist() * pbpair.first;
                    _p += (_lambda_prev * _v);
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;
            }
        }

    private :

        double _L;
        VT _p;
        VT _p0;
        VT _v;
        NT _lambda_prev;
        MT _AA;
        update_parameters _update_parameters;
        VT _lambdas;
        VT _Av;
    };

};


#endif



