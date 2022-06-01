/// This class handles correlation matrices
/// @tparam NT Numeric Type

template<typename Point>
class correl_matrix {
    public:

    /// The numeric/matrix/vector types we use
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
   
    /// The linear matrix inequality that describes the spectrahedron
    LMI<NT, MT, VT> lmi;

    double maxDouble = std::numeric_limits<double>::max();

    /// The correlation matrix
    MT matrix;
    
    /// The coefficients of the correlation matrix
    VT coeffs;
    
    /// The size of the matrix
    unsigned int n;

    /// The dimension of the vector x
    unsigned int d;

    EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;

    VT grad;
    std::pair<PointType, NT> _inner_ball;


    Correlation_matrix(int n, VT coeffs){
        this->n = n;
        d = n*(n-1)/2;
        this->coeffs = coeffs;
        int i, j;
        for(i = 0; i < n; ++i){
            for(j = 0; j < i; ++j){
                matrix(i,j) = coeffs[];
            }
            matrix(i,i) = 1;
        }
        for(i = 0; i < n; ++i){
            for(j = i+1; j < n; ++j){
                matrix(j,i) = matrix(i,j);
            }
        }
    }

    // A small function to build a correlation matrix from a vector of entries
    MT rebuildMatrix(const Point &p, const unsigned int n){
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

    /// \returns The dimension of vector x
    unsigned int dimension() const {
        return d;
    }

    /// \return The matrices A0, A1, ..., Ad
    MT getMatrix() const {
        return matrix;
    }

    /// \returns The size of the matrix
    unsigned int getSizeOfMatrix() const {
        return n;
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    MT* const getMatrix(const int i) {
        return &matrix;
    }


    /// Prints the correlation matrix
    void print() const {
        std::cout << matrix << std::endl;
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

    std::pair<PointType, NT> ComputeInnerBall() {
        return std::pair<PointType, NT>(_inner_ball.first, 1/std::sqrt(d));
    }
    // std::pair<PointType, NT> ComputeInnerBall() {

    //     NT radius = maxDouble;
    //     interior_point = _inner_ball.first.getCoefficients();
    //     for (unsigned int i = 0; i < dimension(); ++i) {

    //         std::pair<NT, NT> min_max = EigenvaluesProblem.symGeneralizedProblem(precomputedValues.A, *(lmi.getMatrix(coordinate)));

    //         if (min_max.first < radius) radius = min_max.first;
    //         if (-min_max.second < radius) radius = -min_max.second;
    //     }

    //     radius = radius / std::sqrt(NT(dimension()));
    //     _inner_ball.second = radius;

    //     return std::pair<PointType, NT>(_inner_ball.first, radius);
    // }
};