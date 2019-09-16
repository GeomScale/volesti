// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis


#ifndef LOW_DIM_SAMPLING_H
#define LOW_DIM_SAMPLING_H

template <class MT, class VT>
std::pair<MT,VT> get_basis(MT Aeq, VT beq){

    MT Ainv = A.completeOrthogonalDecomposition().pseudoInverse();
    int d = A.cols();

    MT AA = MT::Identity(d,d) - Ainv * A;
    Eigen::ColPivHouseholderQR<MT> qrdecomp(AA);
    MT Q = qrdecomp.householderQ();

    Eigen::FullPivLU<MT> lu_decomp(SS);
    int rank = lu_decomp.rank();
    MT retMat(d,rank);

    for (int i = 0; i < rank; ++i) {
        retMat.col(i) = AA.col(i);
    }

    return std::pair<MT,VT> (retMat, Ainv * beq);

}

template <class Hpolytope, class MT, class VT>
Hpolytope get_poly_and_mat_transform(MT A, VT b, MT Aeq, VT beq){

    std::pair<MT,VT> basis = get_basis(Aeq, beq);

    


}


#endif