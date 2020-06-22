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


#ifndef GET_FULL_DIMENSIONAL_POLYTOPE
#define GET_FULL_DIMENSIONAL_POLYTOPE

template <typename Hpolytope, typename Polytope>
std::pair<Hpolytope, std::pair<MT, VT> > get_full_dimensional_polytope(Polytope &P)
{
    typedef typename Hpolytope::NT NT;

    MT Aeq = P.get_Aeq();
    MT A = P.get_mat();
    VT b = P.get_vec();
    VT b = P.get_beq();

    FullPivLU<MatrixXd> lu(Aeq);
    MatrixXd N = lu.kernel();

    VT p = Aeq.colPivHouseholderQr().solve(beq);
    b = b - A * p;
    A = A * N;

    Hpolytope HP(P.dimension(), A, b);

    return std::pair<Hpolytope, std::pair<MT, VT> >(HP, std::pair<N, p>);

}

#endif
