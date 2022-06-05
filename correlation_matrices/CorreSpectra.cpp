/// This class handles the spectrahedra of correlation matrices
/// @tparam Point
template<typename Point>
class CorreSpectra {
    public:

    /// The numeric/matrix/vector types
    typedef typename Point::FT                                NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
   
    /// The size of the matrix
    unsigned int n;

    /// The dimension of the vector x
    unsigned int d;

    /// The linear matrix inequality that describes the spectrahedron
    std::vector<MT> lmi;

    std::pair<Point, NT> _inner_ball;

    /// The gradient vector
    VT grad;

    EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;

    CorreSpectra(unsigned int n){
        this->n = n;
        d = n*(n-1)/2;
        MT A;
        lmi.push_back(MT::Identity(n, n));
        for(i = 0; i < n; ++i){
            for(j = i+1; j < n; ++j){
                A = MT::Zero(n, n);
                A(i,j) = 1;
                A(j,i) = 1;
                list_Mat.push_back(A);
            }
        }
        _inner_ball.first = Point(d);
        _inner_ball.second = 1/std::sqrt(d);
    }

    /// \returns The dimension of vector x
    unsigned int dimension() const {
        return d;
    }

    /// \returns The size of the matrix
    unsigned int getSizeOfMatrix() const {
        return n;
    }
    
    /// \return The matrices A0, A1, ..., Ad
    std::vector<MT> getMatrices() const {
        return lmi;
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    MT* const getMatrix(const int i) {
        return &lmi.at(i);
    }

    std::pair<int,int> getMatrixIndices(int n, int ind){
        int i = 2*n-1-sqrt((2*n-1)*(2*n-1) - 8*ind)/2;
        int j = ind - (2*n-1-i)*i/2;
        return new pair(i,j);
    }

    /// Build a correlation matrix from a vector of entries
    MT buildMatrix(const Point &p, const unsigned int n){
        VT xvector = p.getCoefficients();
        MT A = MT::Identity(n,n);
        double coeff;
        for(int i = 0; i < n ; ++i){
            for(int j = i+1; j < n; ++j){
                int ind = ((((n<<1)-i-2)*(i+1)) >> 1)  + j - n;
                coeff = xvector[ind];
                A(i,j) = coeff;
                A(j,i) = coeff;
            }
        }
        return A;
    }

    /// Find out is lmi(current position) = mat is in the exterior of the spectrahedron
    /// \param mat a matrix where mat = lmi(current position)
    /// \return true if position is outside the spectrahedron
    
    bool is_in(Point const& p, NT tol=NT(0))
    {   
        EigenvaluesProblems<NT, MT, VT> eigs;
        NT eival = eigs.findSymEigenvalue(matrix);
        return eival > 0;
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(VT p, VT const& e, VT &ret) const {
        NT sum_sqqrt_sq = NT(0);
        for (int ind = 0; ind < d; ++ind) {
            std::pair<int,int> indices = getMatrixIndices(n, ind);
            ret(i) = e[indices.first]*e[indices.second];
            sum_sqqrt_sq += ret(ind)*ret(ind);
        }
        ret /= std::sqrt(sum_sqqrt_sq); //normalize
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] r A point on the boundary of the spectrahedron
    /// \param[in] v The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    template <typename update_parameters>
    void compute_reflection(Point &v, Point const& r, update_parameters& ) const 
    {
        VT grad(d);
        lmi.normalizedDeterminantGradient(r.getCoefficients(), precomputedValues.eigenvector, grad);

        // v: original direction s: the surface normal
        // reflected direction = v - 2 <v,s>*s
        NT dot = 2 * v.dot(grad);
        v += -dot * Point(grad);
    }

    std::pair<double, int> intersection(spectrahedron &P, const Point &x, const Point &v, const unsigned int k){
        double tau, tmp;
        int j = 0;
        if(v[0] > 0){
            tau = (1-x[0])/v[0];   
        }else{
            tau = -(1 + x[0])/v[0];
        }
        for(int i = 1; i < k; ++i){
            if(v[i] > 0){
                tmp = (1 - x[i])/v[i];
            }else{
                tmp = -(1 + x[i])/v[i];
            }
            if(tau > tmp){
                tau = tmp;
                j = i;
            }
        }
        tmp = P.positiveLinearIntersection(x.getCoefficients(), v.getCoefficients());
        if(tau > tmp){
            tau = tmp;
            j = -1;
        }
        std::pair<double, int> res(tau,j);
        return res;
    }

    /// Computes the intersection of the line a + tv,
    /// assuming we start at t=0 and that b has zero everywhere and 1 in its i-th coordinate.
    /// Solve the generalized eigenvalue problem A+tB = lmi(a + tv)
    /// So A = lmi(a) and B=lmi(v) - A0
    /// When 
    /// \param[in] a Input vector
    /// \param[in] coordinate Indicator of the i-th coordinate, 1 <= coordinate <= dimension
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT coordinateIntersection(int const coordinate) {

/// In this case, it is very special
        return EigenvaluesProblem.symGeneralizedProblem(matrix, *.getMatrix(coordinate)));
    }

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
};