/// Find the minimum positive and maximum negative eigenvalues of the generalized eigenvalue
    /// problem lB - A, where A, B symmetric and A positive definite.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \return The pair (minimum positive, maximum negative) of eigenvalues
    NTpair symGeneralizedProblem(MT const & A, MT const & B) {

        int matrixDim = A.rows();

        // Spectra solves Xv=lYv, where Y positive definite
        // Set X = B, Y=-A. Then, the eigenvalues we want are the minimum negative
        // and maximum positive eigenvalues of Xv=lYv.

        // Construct matrix operation object using the wrapper classes provided by Spectra
        Spectra::SparseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(-A);

        // Construct generalized eigen solver object
        // requesting the minmum negative and largest positive eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::BOTH_ENDS, Spectra::SparseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 5 < matrixDim ? 5 : matrixDim);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        if (geigs.info() != Spectra::SUCCESSFUL)
            return {NT(0), NT(0)};

        Eigen::VectorXd evalues;
        double lambdaMinPositive, lambdaMaxNegative;

        evalues = geigs.eigenvalues();

        // get the eigenvalues of the original problem
        lambdaMinPositive = 1 / evalues(0);
        lambdaMaxNegative = 1 / evalues(1);

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    NT minPosLinearEigenvalue(MT const & A, MT const & B, VT &eigvec) 
    {
        int matrixDim = A.rows();
        double lambdaMinPositive;

        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(-A);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::LARGEST_ALGE,  Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY> 
            geigs(&op, &Bop, 1, 15 < matrixDim ? 15 : matrixDim);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        VT evalues;
        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
            eigvec = geigs.eigenvectors().col(0);
        }

        lambdaMinPositive = 1 / evalues(0);

        return lambdaMinPositive;
    }
    /// Find the smallest eigenvalue of M
    /// \param M a symmetric matrix
    /// \return smallest eigenvalue
    NT findSymEigenvalue(MT const & M) {
        EigenDenseMatrix<NT> _M(&M);

        // This parameter is for Spectra. It must be larger than #(requested eigenvalues) + 2
        // and smaller than the size of matrix;
        int ncv = M.cols()/10 + 5;
        if (ncv > M.cols()) ncv = M.cols();

        Spectra::SymEigsSolver<NT, Spectra::LARGEST_ALGE, EigenDenseMatrix<NT> > eigs(&_M, 1, ncv);
        // compute
        eigs.init();
        eigs.compute(50000);
        if(eigs.info() == Spectra::SUCCESSFUL) {
            return eigs.eigenvalues()(0);
        }
        else {
            std::cout << "Spectra failed\n";
            return NT(0);
        }
    }