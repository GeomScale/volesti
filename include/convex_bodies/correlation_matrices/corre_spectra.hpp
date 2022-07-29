#ifndef VOLESTI_CORRE_SPECTRAHEDRON_H
#define VOLESTI_CORRE_SPECTRAHEDRON_H

template <typename NT, typename MT, typename VT>
struct Precompute {

    /// These flags indicate whether the corresponding matrices are computed
    bool computed_A = false;
    bool computed_B = false;

    /// The matrices the method positiveIntersection receives from its previous call
    /// if the flag first_positive_intersection is true.
    /// Matrix A is also used in coordinateIntersection
    MT A, B;

    /// In method positive_intersect, the distance we are computing corresponds
    /// to the minimum positive eigenvalue of a quadratic eigenvalue problem.
    /// This will hold the eigenvector for that eigenvalue
    VT eigenvector;

    /// Sets all flags to false
    void resetFlags() {
        computed_A = computed_B = false;
    }

    void set_mat_size(int const& n) 
    {
        A = -MT::Identity(n,n);
        B.setZero(n, n);
        eigenvector.setZero(n);
    }
};

/// This class handles the spectrahedra of correlation matrices
/// @tparam Point
template<typename Point>
class CorreSpectra : public Spectrahedron<Point> {
    public:

    /// The numeric/matrix/vector types we use
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
    typedef Precompute<NT, MT, VT> PrecomputationOfValues;

    /// The size of the matrix
    unsigned int n;

    PrecomputationOfValues _precomputedValues;

    /// Constructor of correlation matrix spectrahedra
    /// \param[in] : matrix size
    CorreSpectra(unsigned int n){
        int i,j;
        this->n = n;
        this->d = n*(n-1)/2;
        this->_inner_ball.first = PointType(this->d);
        this->_inner_ball.second = 1/std::sqrt(this->d);
        _precomputedValues.set_mat_size(n);
    }

    /// \returns The size of the matrix
    unsigned int matrixSize() const {
        return n;
    }

    std::pair<PointType, NT> getInnerBall() const {
        return this->_inner_ball;
    }

    /// Build a correlation matrix from a vector of entries
    /// \param[in] vector of coefficients
    /// \param[in] the matrix to be assigned
    void buildMatrix(const VT &pvector, MT & mat){
        NT coeff;
        int i, j, ind = 0;
        for(i = 0; i < n ; ++i){
            mat(i,i) = -1;
        }
        for(i = 0; i < n ; ++i){
            for(j = i+1; j < n; ++j){
                coeff = -pvector[ind];
                mat(i,j) = mat(j,i) = coeff;
                ++ind;
            }
        }
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] r A point on the boundary of the spectrahedron
    /// \param[in] v The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    template <typename update_parameters>
    void compute_reflection(PointType &v, PointType const& r, update_parameters& ) const {   
        VT grad(this->d);
        VT e = _precomputedValues.eigenvector;
        int i, j, ind = 0;
        NT sum_sq = NT(0);

        for(i = 0; i < n ; ++i){
            for(j = i+1; j < n; ++j){
                grad(ind) = e[i]*e[j];
                sum_sq += grad(ind)*grad(ind);
                ++ind;
            }
        }
        NT dot = v.dot(grad);
        dot = 2 * dot / sum_sq;
        v -= dot * PointType(grad);
    }

    /// Construct the generalized eigenvalue problem \[Bt - A \] for positive_intersect.
    /// \param[in] p Input vector
    /// \param[in] v Input vector
    /// \param[in, out] _precomputedValues Holds matrices B = I - A(v), A = A(p)
    void createMatricesForPositiveLinearIntersection(const VT& p, const VT& v) {
        if (true) {
            VT pvector = p, vvector = v;
            NT coeff;
            int i, j, ind =0;
            for(i = 0; i < n ; ++i){
                for(j = i+1; j < n; ++j){
                    coeff = -pvector[ind];
                    _precomputedValues.A(i,j) = _precomputedValues.A(j,i) = coeff;
                    coeff = -vvector[ind];
                    _precomputedValues.B(i,j) = _precomputedValues.B(j,i) = coeff;
                    ++ind;
                }
            }
            _precomputedValues.computed_B = true;
        }
    }

    NT positiveLinearIntersection(VT const & p, VT const & v) {
        createMatricesForPositiveLinearIntersection(p, v);
        return this->EigenvaluesProblem.minPosLinearEigenvalue2(-_precomputedValues.A, _precomputedValues.B,
                                                                _precomputedValues.eigenvector);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v) {
        NT pos_inter = positiveLinearIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT& ,
                                               NT const&){
        return line_positive_intersect(r, v);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT&){
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT& ,
                                               NT const&,
                                               update_parameters&)
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT&,
                                               NT const&,
                                               MT const&,
                                               update_parameters& )
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType const& r,
                                                     PointType const& v,
                                                     VT&,
                                                     VT&,
                                                     update_parameters&)
    {
        return line_positive_intersect(r, v);
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: the eigenvector A(p)*e = 0
    /// \param[out] ret The unit normal vector at p
    void unit_normal(VT p, VT const& e, VT &ret) const {
        int i, j, ind = 0;
        NT sum_sqqrt_sq = NT(0);
        for(i = 0; i < n ; ++i){
            for(j = i+1; j < n; ++j){
                ret(ind) = e[i]*e[j];
                sum_sqqrt_sq += ret(ind)*ret(ind);
                ++ind;
            }
        }
        ret /= std::sqrt(sum_sqqrt_sq); //normalize
    }

    /// Test if a point p is in the spectrahedron
    /// \param p is the current point
    /// \return true if position is outside the spectrahedron
    int is_in(PointType const& p, NT tol=NT(0)) {  
        return !isExterior(p.getCoefficients());
    }

    bool isExterior(VT const & pos) {
        if(!_precomputedValues.computed_A){
            buildMatrix(pos, n, _precomputedValues.A);
        }
        return isExterior(_precomputedValues.A);
    }

    bool isExterior(MT const & mat) {
        return !this->EigenvaluesProblem.isPositiveSemidefinite(-mat);
    }
};

#endif //VOLESTI_CORRE_SPECTRAHEDRON_H