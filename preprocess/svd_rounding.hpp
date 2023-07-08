// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef SVD_ROUNDING_HPP
#define SVD_ROUNDING_HPP


template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename Point,
    typename MT,
    typename VT,
    typename RandomNumberGenerator
>
void svd_on_sample(Polytope &P, Point &p, unsigned int const& num_rounding_steps, MT &V, VT &s, VT &Means,
                   unsigned int const& walk_length, RandomNumberGenerator &rng)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    typedef RandomPointGenerator <walk> RandomPointGenerator;
    PushBackWalkPolicy push_back_policy;

    unsigned int N = num_rounding_steps;

    std::list<Point> randPoints;
    MT RetMat(N, P.dimension());
    RandomPointGenerator::apply(P, p, N, walk_length, randPoints,
                                push_back_policy, rng);

    int jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        RetMat.row(jj) = (*rpit).getCoefficients().transpose();
    }

    for (int i = 0; i < P.dimension(); ++i) {
        Means(i) = RetMat.col(i).mean();
    }

    for (int i = 0; i < N; ++i) {
        RetMat.row(i) = RetMat.row(i) - Means.transpose();
    }

    Eigen::BDCSVD<MT> svd(RetMat, Eigen::ComputeFullV);
    s = svd.singularValues() / svd.singularValues().minCoeff();

    if (s.maxCoeff() >= 2.0) {
        for (int i = 0; i < s.size(); ++i) {
            if (s(i) < 2.0) {
                s(i) = 1.0;
            }
        }
        V = svd.matrixV();
    } else {
        s = VT::Ones(P.dimension());
        V = MT::Identity(P.dimension(), P.dimension());
    }
}


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

std::tuple<MT, VT, NT> svd_rounding(Polytope &P,
                                    std::pair<Point,NT> &InnerBall,
                                    const unsigned int &walk_length,
                                    RandomNumberGenerator &rng)
{
    NT tol = 0.00000001;
    NT R = std::pow(10,10), r = InnerBall.second;

    int n = P.dimension(), m = P.num_of_hyperplanes();

    Polytope P_old(P);

    MT old_A(m,n), A = P.get_mat(), T = MT::Identity(n,n), round_mat, r_inv;
    VT T_shift = VT::Zero(n), shift(n), s(n);

    Point p(n);

    bool done = false, last_round_under_p, fail;

    unsigned int tries=0, num_rounding_steps = 10 * n, rounding_samples = 0, round_it;
    NT max_s, s_cutof, p_cutof, num_its, prev_max_s = std::numeric_limits<NT>::max(),
       s_cutoff, p_cutoff;
    MT V(n,n), S(n,n);

    while (!done) {

        T = MT::Identity(n, n);
        T_shift = VT::Zero(n);

        round_it = 1;
        max_s = std::numeric_limits<NT>::max();
        s_cutoff = 2.3;
        p_cutoff = 10.0;
        last_round_under_p = false;
        fail = false;
        num_its = 20;

        while (max_s > s_cutoff && round_it <= num_its) {

            p = InnerBall.first;
            svd_on_sample<WalkTypePolicy>(P, p, num_rounding_steps, V, s,
                                          shift, walk_length, rng);

            rounding_samples = rounding_samples + num_rounding_steps;
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
                break;
            }

            round_it++;
            prev_max_s = max_s;

            P.shift(shift);
            P.linear_transformIt(round_mat);
            InnerBall = P.ComputeInnerBall();
            T_shift += T * shift;
            T = T * round_mat;
        }
        if (round_it <= num_its && !fail) {
            done = true;
        } else {
            tries = tries + 1;
            num_rounding_steps = num_rounding_steps * 2.0;
            P = P_old;
        }
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(T.determinant()));
    return result;
}


#endif
