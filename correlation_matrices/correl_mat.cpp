/// This class handles correlation matrices
/// @tparam NT Numeric Type

template<typename Point>
class CorreSpectra {
    public:

    /// The numeric/matrix/vector types we use
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
   
    /// The size of the matrix
    unsigned int n;

    /// The dimension of the vector x
    unsigned int d;

    /// The linear matrix inequality that describes the spectrahedron
    std::vector<MT> lmi;

    std::pair<PointType, NT> _inner_ball;

    EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;

    VT grad;

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
        return matrix;
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    MT* const getMatrix(const int i) {
        return &lmi[i];
    }

    /// Find out is lmi(current position) = mat is in the exterior of the spectrahedron
    /// \param mat a matrix where mat = lmi(current position)
    /// \return true if position is outside the spectrahedron
    
    bool is_in(PointType const& p, NT tol=NT(0))
    {   
        EigenvaluesProblems<NT, MT, VT> eigs;
        NT eival = eigs.findSymEigenvalue(matrix);
        return eival > 0;
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

    void reflection(spectrahedron P, Point &p, Point &v, const int flag){
        if(flag != -1){
            v.set_coord(flag, - v.getCoefficients()(flag));
            return;
        }
        P.compute_reflection(v, p, flag);
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(VT r, VT const& e, VT &ret) const {
        NT* ret_data = ret.data();
        NT sum_sqqrt_sq = NT(0);
        for (int i = 0; i < d; i++) {
            // todo, use iterators
            *ret_data = e.dot(matrices[i+1].template selfadjointView< Eigen::Lower >() * e);
            sum_sqqrt_sq += (*ret_data) * (*ret_data);
            ret_data++;
        }

        //normalize
        ret /= std::sqrt(sum_sqqrt_sq);
    }

    /// Computes the distance d one must travel on the line a + tb,
    /// assuming we start at t=0 and that b has zero everywhere and 1 in its i-th coordinate.
    /// We must solve the generalized eigenvalue problem A+tB, where A = lmi(a) and B=(lmi) - A0 = A_i
    /// If the flag precomputedValues,computed_A is true, the matrix A is not computed.
    /// \param[in] a Input vector
    /// \param[in] coordinate Indicator of the i-th coordinate, 1 <= coordinate <= dimension
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT coordinateIntersection(int const coordinate) {
        return EigenvaluesProblem.symGeneralizedProblem(matrix, *.getMatrix(coordinate)));
    }
};