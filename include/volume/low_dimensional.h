// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef LOW_DIMENSIONAL_H
#define LOW_DIMENSIONAL_H

template <class Polytope, class MT, class VT>
std::pair<Polytope, VT> get_low_dimensional_poly(MT A, VT b, MT Aeq, VT beq, MT &W) {

    typedef typename Polytope::NT 	NT;
    Polytope HP;
    int m = Aeq.rows(), n = Aeq.cols();

    VT x = Aeq.fullPivLu().solve(beq);
    b = b - A*x;

    Eigen::JacobiSVD<MT> svd(Aeq.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    MT T2 = svd.matrixU().transpose();

    W.resize(n,n-m);
    for (int i = 0; i < n-m; ++i) {
        W.col(i) = T2.col(m+i);
    }

    MT A2 = A*W;
    HP.init(n-m, A2, b);
    return std::pair<Polytope, VT>(HP,x);

}

#endif
