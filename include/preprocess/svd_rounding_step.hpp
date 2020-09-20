// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef SVD_ROUNDING_SINGLE_STEP_HPP
#define SVD_ROUNDING_SINGLE_STEP_HPP


template
<
    typename WalkTypePolicy,
    typename MT,
    typename VT,
    typename Polytope,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>

void svd_rounding_single_step(Polytope &P,
                                    std::pair<Point,NT> &InnerBall,
                                    const unsigned int &walk_length,
                                    MT &T, VT &T_shift,
                                    unsigned int &num_rounding_steps,
                                    bool &fail, bool &converged, bool &last_round_under_p,
                                    NT &max_s, unsigned int &round_it
                                    RandomNumberGenerator &rng)
{
    //NT tol = 0.00000001;
    //NT R = std::pow(10,10), r = InnerBall.second;

    int n = P.dimension();//, m = P.num_of_hyperplanes();

    //Polytope P_old(P);

    MT round_mat, r_inv;
    VT shift(n), s(n);

    Point p(n);

    //bool last_round_under_p;

    //unsigned int round_it;
    NT num_its, prev_max_s = std::numeric_limits<NT>::max(),
       s_cutoff, p_cutoff;
    MT V(n,n), S(n,n);

    round_it = 1;
    max_s = std::numeric_limits<NT>::max();
    s_cutoff = 2.3;
    p_cutoff = 10.0;
    //last_round_under_p = false;
    num_its = 20;

    //while (max_s > s_cutoff && round_it <= num_its) {

        p = InnerBall.first;
        svd_on_sample<WalkTypePolicy>(P, p, num_rounding_steps, V, s,
                                      shift, walk_length, rng);

        //rounding_samples = rounding_samples + num_rounding_steps;
        max_s = s.maxCoeff();

        if (max_s <= p_cutoff && max_s > s_cutoff) {
            if (last_round_under_p) {
                num_rounding_steps = num_rounding_steps * 2;
                p = InnerBall.first;
                svd_on_sample<WalkTypePolicy>(P, p, num_rounding_steps, V, s,
                                              shift, walk_length, rng);
                max_s = s.maxCoeff();
            } else {
                last_round_under_p = true;
            }
        } else {
            last_round_under_p = false;
        }
        S = s.asDiagonal();
        round_mat = V * S;
        r_inv = VT::Ones(n).cwiseProduct(s.cwiseInverse()).asDiagonal() * V.transpose();

        if (round_it != 1 && max_s >= NT(4) * prev_max_s) {
            fail = true;
            //break;
        }

        round_it++;
        prev_max_s = max_s;

        P.shift(shift);
        P.linear_transformIt(round_mat);
        //InnerBall = P.ComputeInnerBall();
        T_shift += T * shift;
        T = T * round_mat;

        if (max_s <= s_cutoff || round_it > num_its) {
            converged = true;
        }
    //}
        

    //std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(T.determinant()));
    //return result;
}


#endif

