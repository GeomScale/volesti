// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_NONCONVEX_GCW_WALK_HPP
#define RANDOM_WALKS_UNIFORM_NONCONVEX_GCW_WALK_HPP


#include "sampling/sphere.hpp"
#include <cmath>

// Random directions hit-and-run walk with uniform target distribution

struct GCWalkOpt
{
    struct parameters {};
    parameters param;

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::NT NT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p, RandomNumberGenerator& rng)
    {
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p,
         RandomNumberGenerator& rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      VT& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            GetDirectionTangentPlane<VT>::apply(p, _v, rng);
            //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            _bpair = P.gc_intersect_all_roots(p, _v, _lamdas, _Av,
                                                       _lambda);
            _lambda = pick_next_angle(_bpair.first, _bpair.second, rng, P, p, _v);
            p = (cos(_lambda) * p) + (sin(_lambda) * _v);
            VT q = P.get_mat()*p - P.get_vec();
            for (int i=0; i<P.num_of_hyperplanes(); i++)
            {
                if (q(i)>NT(0))
               {
                    std::cout<<"outside from sampling, q: "<<q(i)<<std::endl;
                    exit(-1);
                }
            }
        }
        //p = _p;
    }



    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           VT& p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _v.setZero(P.dimension());

        std::cout<<"Initialize..."<<std::endl;
        //std::cout<<"_v = "<<_v.transpose()<<"\n"<<std::endl;
        if (P.is_in(p) == 0) {
                std::cout<<"initial point out"<<std::endl;
                //std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<", index_neg = "<<index_neg<<std::endl;
                exit(-1);
            }

        GetDirectionTangentPlane<VT>::apply(p, _v, rng);
        std::cout<<"_v = "<<_v.transpose()<<"\n"<<std::endl;
        std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
        _bpair = P.gc_intersect_all_roots(p, _v, _lamdas, _Av);
        _lambda = pick_next_angle(_bpair.first, _bpair.second, rng, P, p, _v);
        p = (cos(_lambda) * p) + (sin(_lambda) * _v);
        VT q = P.get_mat()*p - P.get_vec();
        for (int i=0; i<P.num_of_hyperplanes(); i++)
        {
            if (q(i)>NT(0))
            {
                std::cout<<"outside from sampling, q: "<<q(i)<<std::endl;
                exit(-1);
            }
        }
    }

    template <typename BallPolytope>
    inline NT pick_next_angle(VT const& neg_roots, VT const& pos_roots, RandomNumberGenerator &rng, BallPolytope const& P, VT const& p, VT const& v) {
        
        int len_pos = pos_roots.rows(), len_neg = neg_roots.rows();

        int index_pos = 0, index_neg = 0;
        std::vector< std::pair<NT, NT> > segments;
        std::vector<NT> lengths;   
        NT theta1, theta2, theta_temp, sum_lenghts = NT(0);
        std::cout<<"len_pos = "<<len_pos<<std::endl;
        std::cout<<"len_neg = "<<len_neg<<std::endl;
        VT q;

        if (len_pos == 1 && len_neg == 1) 
        {
            return rng.sample_urdist() * (pos_roots[0] - neg_roots[0]) + neg_roots[0];
        } 
        else
        {
            theta2 = pos_roots[0];
            theta1 = neg_roots[len_neg - 1];
            segments.push_back(std::pair<NT, NT>(theta1, theta2));
            lengths.push_back(theta2 - theta1);
            sum_lenghts += (theta2 - theta1);
            theta_temp = (theta2 + theta1) / NT(2);
            q = (cos(theta_temp) * p) + (sin(theta_temp) * v);
            if (P.is_in(q) == 0) {
                std::cout<<"point out"<<std::endl;
                std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<std::endl;
                exit(-1);
            }
        }
        std::cout<<"iterating over positives roots"<<std::endl;
        while (index_pos < len_pos)
        {
            if (index_pos == (len_pos - 1)) {
                std::cout<<"index_pos = "<<index_pos<<std::endl;
                index_neg = 1;
                theta_temp = neg_roots[0] + NT(2)*M_PI;
                theta1 = std::min(theta_temp, pos_roots[index_pos]);
                theta2 = std::max(theta_temp, pos_roots[index_pos]);
                std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<std::endl;
            } else {
                theta1 = pos_roots[index_pos];
                theta2 = pos_roots[index_pos + 1];
            }
            segments.push_back(std::pair<NT, NT>(theta1, theta2));
            lengths.push_back(theta2 - theta1);
            sum_lenghts += (theta2 - theta1);
            
            theta_temp = (theta2 + theta1) / NT(2);
            q = (cos(theta_temp) * p) + (sin(theta_temp) * v);
            if (P.is_in(q) == 0) {
                std::cout<<"point out"<<std::endl;
                std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<", index_pos = "<<index_pos<<std::endl;
                if (index_pos < (len_pos-1)){
                    theta_temp = (pos_roots[index_pos+1] + pos_roots[index_pos+2]) / NT(2);
                    q = (cos(theta_temp) * p) + (sin(theta_temp) * v);
                    if (P.is_in(q) == 0) {
                        std::cout<<"also next seg point out"<<std::endl;
                    }
                }
                //exit(-1);
            }
            index_pos += 2;
            std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<std::endl;
        }
        std::cout<<"iterating over negatives roots"<<std::endl;
        while (index_neg < len_neg)
        {
            if (index_neg == (len_neg - 1)){
                break;
            }
            theta1 = neg_roots[index_neg];
            theta2 = neg_roots[index_neg + 1];
            segments.push_back(std::pair<NT, NT>(theta1, theta2));
            lengths.push_back(theta2 - theta1);
            sum_lenghts += (theta2 - theta1);
            
            theta_temp = (theta2 + theta1) / NT(2);
            q = (cos(theta_temp) * p) + (sin(theta_temp) * v);
            if (P.is_in(q) == 0) {
                std::cout<<"point out"<<std::endl;
                std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<", index_neg = "<<index_neg<<std::endl;
                //exit(-1);
                if (index_neg < (len_neg-1)){
                    theta_temp = (neg_roots[index_neg+1] + neg_roots[index_neg+2]) / NT(2);
                    q = (cos(theta_temp) * p) + (sin(theta_temp) * v);
                    if (P.is_in(q) == 0) {
                        std::cout<<"also next seg point out"<<std::endl;
                    }
                }
            }
            index_neg += 2;
            std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<std::endl;
            //std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<", theta2 - theta1 = "<<theta2-theta1<<std::endl;
        }

        std::cout<<"num_segs = "<<segments.size()<<", num_lengths = "<<lengths.size()<<std::endl;
        NT sum_lens2 = 0;
        for (int i=0; i<lengths.size(); i++)
        {
            lengths[i] = lengths[i] / sum_lenghts;
            std::cout<<lengths[i]<<", ";
            sum_lens2 += lengths[i];
            
        }
        std::cout<<"\n";
        std::cout<<"sum_lens2 = "<<sum_lens2<<std::endl;

        return get_next_theta(segments, lengths, rng);
    }


    inline NT get_next_theta(std::vector< std::pair<NT, NT> > &segments, std::vector<NT> &lengths, RandomNumberGenerator &rng)
    {
        NT u = rng.sample_urdist(), sum_lens = NT(0);
        int n_segs = segments.size(), index;

        for(int i=0; i<n_segs; i++)
        {
            sum_lens += lengths[i];
            if (u < sum_lens)
            {
                index = i;
                break;
            }
        }

        NT theta1 = segments[index].first, theta2 = segments[index].second;
        std::cout<<"theta1 = "<<theta1<<", theta2 = "<<theta2<<std::endl;
        NT theta = rng.sample_urdist() * (theta1 - theta2) + theta2;
        std::cout<<"theta = "<<theta<<std::endl;

        return theta;
    }

private :

    //Point _p;
    NT _lambda;
    VT _lamdas;
    VT _Av;
    VT _v;
    VT _neg_roots;
    VT _pos_roots;
    std::pair<VT, VT> _bpair; 
};

};


#endif // RANDOM_WALKS_UNIFORM_GCW_WALK_HPP
