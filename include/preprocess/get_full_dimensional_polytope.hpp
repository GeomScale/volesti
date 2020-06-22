
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

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

template <typename H_polytope, typename MT, typename VT>
std::pair<Hpolytope, std::pair<MT, VT> > get_full_dimensional_polytope(MT A, VT b, MT Aeq, VT beq)
{
    typedef typename H_polytope::NT NT;

    VT p = Aeq.colPivHouseholderQr().solve(beq);

    FullPivLU<MT> lu(Aeq);
    MatrixXd N = lu.kernel();

    b = b - A * p;
    A = A * N;

    H_polytope HP;
    HP.init(A.cols(), A, b);

    return std::pair<Hpolytope, std::pair<MT, VT> >(HP, std::pair<MT,VT>(N, p));

}

#endif
