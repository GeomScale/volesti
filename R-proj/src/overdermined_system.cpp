// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>

//' Internal rcpp function to compute the full dimensional polytope when a low dimensional is given
//'
//' @param P A low dimensional convex polytope in H-representation.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the full dimensional polytope, a numerical matrix of the inverse
//'         linear transformation that is applied on the input polytope, the numerical vector - point that the
//'         input polytope is shifted and the product of the singular values of the matrix of the linear map 
//'         applied on the input polytope.
//'
// [[Rcpp::export]]
Rcpp::NumericVector solve_overdetermined_linear_system(Rcpp::NumericVector row_ind, Rcpp::NumericVector col_ind, 
                                                       Rcpp::NumericVector values, Rcpp::NumericVector b, int m, int d)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> VTint;
    typedef Eigen::SparseMatrix<NT> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> DMT;

    Eigen::SparseQR<MT, Eigen::COLAMDOrdering<int> >   solver;

    typedef Eigen::Triplet<NT> T;
    std::vector<T> tripletList;

    VTint rows = Rcpp::as<VTint>(row_ind), cols = Rcpp::as<VTint>(col_ind);
    VT vals = Rcpp::as<VT>(values), beq = Rcpp::as<VT>(b);

    int n = rows.size();

    MT A(m, d);
    for (int i=0; i<n; i++) {
        tripletList.push_back(T(rows(i), cols(i), vals(i)));
    }


    A.setFromTriplets(tripletList.begin(), tripletList.end());
    
    //DMT Aeq = DMT(A);
    
    //VT x = Aeq.fullPivLu().solve(beq);

    solver.analyzePattern(A); 
    // Compute the numerical factorization 
    solver.factorize(A); 
    //Use the factors to solve the linear system 
    VT y = solver.solve(beq); 

    
    

    return Rcpp::wrap(y);
}



// [[Rcpp::export]]
Rcpp::NumericMatrix null_space_SparseQR(Rcpp::NumericVector row_ind, Rcpp::NumericVector col_ind, 
                                                       Rcpp::NumericVector values, int m, int d)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> VTint;
    typedef Eigen::SparseMatrix<NT> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> DMT;

    Eigen::SparseQR<MT, Eigen::COLAMDOrdering<int> >   solver;

    typedef Eigen::Triplet<NT> T;
    std::vector<T> tripletList;

    VTint rows = Rcpp::as<VTint>(row_ind), cols = Rcpp::as<VTint>(col_ind);
    VT vals = Rcpp::as<VT>(values);

    int n = rows.size();

    MT A(m, d);
    std::cout<<"d = "<<d<<", m = "<<m<<std::endl;
    std::cout<<"rows = "<<rows.transpose()<<"\n cols = "<<cols.transpose()<<std::endl;
    for (int i=0; i<n; i++) {
        tripletList.push_back(T(rows(i), cols(i), vals(i)));
    }

    std::cout<<"triplets done"<<std::endl;
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    MT B = A.transpose();
    std::cout<<"A done"<<std::endl;
    B.makeCompressed();
    //std::cout<<"B = "<<DMT(B)<<std::endl;
    
    solver.setPivotThreshold(0.000001);
    //solver.compute(B); 
    // Compute the numerical factorization 
    solver.analyzePattern(B);
    solver.factorize(B);

    std::cout<<"factorization done"<<std::endl;
    int rnk = solver.rank();
    rnk=502;
    std::cout<<"rank = "<<rnk<<std::endl;

    //MT Q = solver.matrixQ();

    DMT N = DMT(solver.matrixQ()).block(0, rnk-1, d, d-rnk);

    

    return Rcpp::wrap(N);
}



