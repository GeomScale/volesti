//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef ROUND_ISOTROPY_HPP
#define ROUND_ISOTROPY_HPP


template <typename WalkTypePolicy,
          typename Polytope,
          typename Point,
          typename MT, 
          typename VT, 
          typename RandomNumberGenerator>
void round_svd(Polytope &P, Point &p, unsigned int const& num_rounding_steps, MT &V, VT &s, VT &Means,
            unsigned int const& walk_length, RandomNumberGenerator &rng)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;
    
    typedef RandomPointGenerator <walk> RandomPointGenerator;
    PushBackWalkPolicy push_back_policy;

    unsigned int N = num_rounding_steps / (P.dimension()*P.dimension() + 1);

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
    Point q(Means);
    p = p - q;

    P.shift(Means);

    Eigen::JacobiSVD<MT> svd(RetMat, Eigen::ComputeThinU | Eigen::ComputeThinV);
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


template <
        typename WalkTypePolicy,
        typename MT,
        typename VT,
        typename Polytope,
        typename Point,
        typename NT,
        typename RandomNumberGenerator
>
std::pair< std::pair<MT, VT>, NT > round_isotropy(Polytope &P, std::pair<Point,NT> &InnerBall,
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

    unsigned int tries=0, num_rounding_steps=8*n*n*n, rounding_samples=0, round_it;
    NT max_s, s_cutof, p_cutof, num_its, prev_max_s = std::numeric_limits<NT>::max(), s_cutoff, p_cutoff;
    MT V(n,n), S(n,n);

    while (!done) {

        T = MT::Identity(n, n);
        T_shift = VT::Zero(n);

        round_it = 1;
        max_s = std::numeric_limits<NT>::max();
        s_cutoff = 4.0;
        p_cutoff = 8.0;
        last_round_under_p = false;
        fail = false;
        num_its = std::max(std::log(R) / std::log(20.0), 2.0);

        while (max_s > s_cutoff && round_it <= num_its) {

            p.set_to_origin();
            round_svd<WalkTypePolicy>(P, p, num_rounding_steps, V, s, shift, walk_length, rng);

            rounding_samples = rounding_samples + num_rounding_steps;
            max_s = s.maxCoeff();

            if (max_s <= p_cutoff && max_s > s_cutoff) {
                if (last_round_under_p) {
                    num_rounding_steps = num_rounding_steps * 2.0;
                    p.set_to_origin();
                    round_svd<WalkTypePolicy>(P, p, num_rounding_steps, V, s, shift, walk_length, rng);
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

            P.linear_transformIt(round_mat);
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
    return std::pair< std::pair<MT, VT>, NT > (std::pair<MT, VT>(T, T_shift), T.determinant());
}


template <
        typename WalkTypePolicy,
        typename MT,
        typename VT,
        typename Polytope,
        typename Point,
        typename NT,
        typename RandomNumberGenerator
>
std::pair< std::pair<MT, VT>, NT > round_isotropy(Polytope &P, std::pair<Point,NT> &InnerBall,
                                                  const unsigned int &walk_length,
                                                  RandomNumberGenerator &rng)
{
    unsigned int d = P.dimension();
    MT N = MT::Identity(d,d);
    VT shift = VT::Zero(d);
    std::pair< std::pair< std::pair<MT, VT>, std::pair<MT, VT> >, NT > result = round_isotropy(P, InnerBall, 
                                                                                    walk_length, rng, N, shift);
    std::pair< std::pair<MT, VT>, NT > res;
    res.first.first = result.first.first.first;
    res.first.second = result.first.first.second;
    res.second = result.second;

    return res;
}

#endif
