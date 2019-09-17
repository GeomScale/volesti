// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis


#ifndef LOW_DIM_SAMPLING_H
#define LOW_DIM_SAMPLING_H

template <class MT, class VT>
std::pair<MT,VT> get_basis(MT Aeq, VT beq) {

    MT Ainv = A.completeOrthogonalDecomposition().pseudoInverse();
    int d = A.cols();

    MT AA = MT::Identity(d, d) - Ainv * A;
    Eigen::ColPivHouseholderQR <MT> qrdecomp(AA);
    MT Q = qrdecomp.householderQ();

    Eigen::FullPivLU <MT> lu_decomp(SS);
    int rank = lu_decomp.rank();
    MT retMat(d, rank);

    for (int i = 0; i < rank; ++i) {
        retMat.col(i) = AA.col(i);
    }

    return std::pair<MT, VT>(retMat, Ainv * beq);

}

template <class MT, class VT>
MT get_transformation(MT basis, VT translation) {

    int d = translation.size();

    basis.conservativeResize(basis.rows(), basis.cols() + 1);
    basis.col(basis.cols() - 1) = VT::Zero(basis.rows());

    basis.conservativeResize(basis.rows() + 1, basis.cols());
    basis.row(basis.rows() - 1) = VT::Zero(basis.cols());
    basis(basis.rows() - 1, basis.cols() - 1) = 1;

    MT elimHom = MT::Zero((basis.rows()-1, basis.rows()));
    //elimHom.col(basis.rows()-1) = VT::Zero(basis.rows() - 1);
    for (int i = 0; i < basis.rows()-1; ++i) {
        elimHom(i,i) = 1;
    }

    MT trans = MT::Identity(d+1, d+1);
    for (int j = 0; j < d; ++j) {
        trans(i, d) = trans(i);
    }

    return (elimHom * trans * basis);

}

template <class MT, class VT>
std::pair<MT,VT> change2_homogeneous(MT A, VT b) {

    int d = A.cols(), m = A.rows();

    MT A2(m, d - 1);
    for (int i = 0; i < d - 1; ++i) {
        A2.col(i) = A.col(i);
    }
    VT b2(m);
    b2 = b - A.col(d - 1);

    return std::pair<MT, VT>(A2, b2);

}

template <class Hpolytope, class MT, class VT>
std::pair<Hpolytope,MT> get_poly_and_mat_transform(MT A, VT b, MT Aeq, VT beq) {

    std::pair <MT, VT> basis = get_basis(Aeq, beq);

    MT T = get_transformation(basis.first, basis.second);

    MT A2 = A * T;

    basis = change2_homogeneous(A2, b);

    Hpolytope HP;
    HP.init(basis.first.cols(), basis.first, basis.second);

    return std::pair<Hpolytope, MT>(HP, T);
}


#endif