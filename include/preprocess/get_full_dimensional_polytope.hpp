// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef GET_FULL_DIMENSIONAL_POLYTOPE
#define GET_FULL_DIMENSIONAL_POLYTOPE


template <typename H_polytope, typename MT, typename VT>
std::pair<H_polytope, std::pair<MT, VT> > get_full_dimensional_polytope(MT A, VT b, MT Aeq, VT beq)
{
    typedef typename H_polytope::NT NT;

    VT p = Aeq.colPivHouseholderQr().solve(beq);
    int rnk = Aeq.rows(), d = A.cols();
    MT N(d, d-rnk);
    MT Aeqtr = Aeq.transpose();

    Eigen::FullPivHouseholderQR<MT> qr2 = Aeqtr.fullPivHouseholderQr();
    MT N2 = qr2.matrixQ();

    int col = 0;
    for (int i = rnk; i<d; i++) {
        N.col(col) = N2.col(i);
        col++;
    }

    b = b - A * p;
    A = A * N;

    H_polytope HP;
    HP.init(d-rnk, A, b);

    return std::pair<H_polytope, std::pair<MT, VT> >(HP, std::pair<MT,VT>(N, p));

}

#endif
