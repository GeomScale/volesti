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


template <typename WalkTypePolicy, typename Polytope, typename Point, typename MT, typename VT, typename RandomNumberGenerator>
void round_svd(Polytope &P, Point &p, unsigned int const& num_rounding_steps, MT &V, VT &s, VT &Means,
            unsigned int const& walk_length, RandomNumberGenerator &rng)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    typedef RandomPointGenerator <walk> RandomPointGenerator;
    PushBackWalkPolicy push_back_policy;

    //p = starting_point;
    unsigned int N = num_rounding_steps / (P.dimension()*P.dimension() + 1);
    std::cout<<"N = "<<N<<"\n"<<std::endl;

    std::cout<<"p = "<<p.getCoefficients()<<"\n"<<std::endl;

    std::cout<<"is_in = "<<P.is_in(p)<<"\n"<<std::endl;

    std::cout<<"P.A = "<<P.get_mat()<<"\n"<<std::endl;
    std::cout<<"P.b = "<<P.get_vec()<<"\n"<<std::endl;

    std::list<Point> randPoints;
    MT RetMat(N, P.dimension());
    //VT Means(P.dimension());

    //randPoints.clear();
    RandomPointGenerator::apply(P, p, N, walk_length, randPoints,
                                push_back_policy, rng);
    std::cout<<"points sampled"<<"\n"<<std::endl;

    int jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        RetMat.row(jj) = (*rpit).getCoefficients().transpose();
    }
    std::cout<<"RetMat = "<<RetMat<<"\n"<<std::endl;

    for (int i = 0; i < P.dimension(); ++i) {
        Means(i) = RetMat.col(i).mean();
    }
    std::cout<<"Means = "<<Means<<"\n"<<std::endl;
    for (int i = 0; i < N; ++i) {
        RetMat.row(i) = RetMat.row(i) - Means.transpose();
    }
    std::cout<<"RetMat = "<<RetMat<<"\n"<<std::endl;
    Point q(Means);
    std::cout<<"q = "<<q.getCoefficients()<<"\n"<<std::endl;
    p = p - q;

    P.shift(Means);
    std::cout<<"P.A = "<<P.get_mat()<<"\n"<<std::endl;
    std::cout<<"P.b = "<<P.get_vec()<<"\n"<<std::endl;

    Eigen::JacobiSVD<MT> svd(RetMat, Eigen::ComputeThinU | Eigen::ComputeThinV);

    s = svd.singularValues() / svd.singularValues().minCoeff();
    std::cout<<"s = "<<s<<"\n"<<std::endl;

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
    std::cout<<"V = "<<V<<"\n"<<std::endl;

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
    //typedef typename Polytope::MT MT;
    //typedef typename Polytope::VT VT;
    //typedef typename Polytope::NT NT;

    NT tol = 0.00000001;
    //std::pair<VT, NT> res =  compute_max_inner_ball(P.get_mat(), P.get_vec(), 150, tol);

    //P.shift(InnerBall.first.getCoefficients());

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
        std::cout<<"num_its = "<<num_its<<"\n"<<std::endl;
        p.set_to_origin();
        std::cout<<"starting iterations"<<"\n"<<std::endl;

        while (max_s > s_cutoff && round_it <= num_its) {

            //template <typename WalkTypePolicy, typename Polytope, typename Point, typename MT, typename VT, typename RandomNumberGenerator>
            //void round_svd(Polytope &P, Point &p, unsigned int const& num_rounding_steps, MT &V, VT &s, VT &Means,
            //               unsigned int &walk_length, RandomNumberGenerator &rng)
            p.set_to_origin();
            round_svd<WalkTypePolicy>(P, p, num_rounding_steps, V, s, shift, walk_length, rng);
            std::cout<<"round_svd computed"<<"\n"<<std::endl;

            rounding_samples = rounding_samples + num_rounding_steps;
            max_s = s.maxCoeff();
            std::cout<<"max_s = "<<max_s<<"\n"<<std::endl;

            if (max_s <= p_cutoff && max_s > s_cutoff) {
                if (last_round_under_p) {
                    std::cout<<"last_round_under_p"<<"\n"<<std::endl;
                    //fprintf('Seem to be close to round. Doubling number of steps, not restarting.\n');
                    num_rounding_steps = num_rounding_steps * 2.0;
                    p.set_to_origin();
                    round_svd<WalkTypePolicy>(P, p, num_rounding_steps, V, s, shift, walk_length, rng);
                    std::cout<<"round_svd computed"<<"\n"<<std::endl;
                    max_s = s.maxCoeff();
                } else {
                    last_round_under_p = true;
                }
            } else {
                last_round_under_p = false;
            }
            S = s.asDiagonal();
            std::cout<<"S = "<<S<<"\n"<<std::endl;
            round_mat = V * S;
            std::cout<<"round_mat = "<<round_mat<<"\n"<<std::endl;
            r_inv = VT::Ones(n).cwiseProduct(s.cwiseInverse()).asDiagonal() * V.transpose();
            std::cout<<"r_inv = "<<r_inv<<"\n"<<std::endl;

            if (round_it != 1 && max_s >= NT(4) * prev_max_s) {

                std::cout<<"fail"<<"\n"<<std::endl;
                fail = true;
                break;
            }

            round_it++;
            prev_max_s = max_s;

            P.linear_transformIt(round_mat);
            T_shift += T * shift;
            T = T * round_mat;
            std::cout<<"linear maps computed"<<"\n"<<std::endl;
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


#endif
