#include "chrono"
#include "matrix_operations/EigenvaluesProblems.h"

template <typename NT, typename MT, typename VT>
struct PrecomputationOfValues {

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
        A.setZero(n, n);
        B.setZero(n, n);

        eigenvector.setZero(m);
    }
};

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
    int n;

    /// The dimension of the vector x
    int d;

    /// The linear matrix inequality that describes the spectrahedron
    std::vector<MT> lmi;

    std::pair<Point, NT> _inner_ball;

    typedef PrecomputationOfValues<NT, MT, VT> _PrecomputationOfValues;

    _PrecomputationOfValues precomputedValues;

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
    int dimension() const {
        return d;
    }

    /// \returns The size of the matrix
    int matrixSize() const {
        return n;
    }
    
    /// \return The matrices A0, A1, ..., Ad
    std::vector<MT> getLMI() const {
        return lmi;
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    MT* getLMI(const int i) const {
        return &lmi.at(i);
    }

    /// Build a correlation matrix from a vector of entries
    void buildMatrix(const VT &pvector, const unsigned int n, MT & mat){
        NT coeff;
        int i, j, ind = 0;
        for(int i = 0; i < n ; ++i){
            for(int j = i+1; j < n; ++j){
                // int ind = ((((n<<1)-i-2)*(i+1)) >> 1)  + j - n;
                coeff = xvector[ind];
                mat(i,j) = mat(j,i) = coeff;
                ++ind;
            }
        }
    }

    std::pair<int,int> getMatrixIndices(int n, int ind){
        int i = 2*n-1-sqrt((2*n-1)*(2*n-1) - 8*ind)/2;
        int j = ind - (2*n-1-i)*i/2;
        return new pair(i,j);
    }

    bool isExterior(MT const & mat) {
        EigenvaluesProblems<NT, MT, VT> eigs;
        NT eival = eigs.findSymEigenvalue(matrix);
        return eival > 0;
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] r A point on the boundary of the spectrahedron
    /// \param[in] v The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    template <typename update_parameters>
    void compute_reflection(Point &v, Point const& r, update_parameters& ) const {   
        VT grad(d);
        unit_normal(r.getCoefficients(), precomputedValues.eigenvector, grad);

        // reflected direction = v - 2 <v,s>*s
        v -= -2 * v.dot(grad) * Point(grad);
    }

    /// Construct the generalized eigenvalue problem \[Bt + A \] for positive_intersect.
    /// \param[in] p Input vector
    /// \param[in] v Input vector
    /// \param[in, out] precomputedValues Holds matrices B = I - A(v), A = A(p)
    void createMatricesForPositiveLinearIntersection(const VT& p, const VT& v) {
        if (!precomputedValues.computed) {
            VT pvector = p.getCoefficients(), vvector = v.getCoefficients();
            precomputedValues.A = MT::Identity(n,n);
            precomputedValues.B = MT::Zero(n,n);
            NT coeff;
            int i, j, ind =0;
            for(i = 0; i < n ; ++i){
                for(j = i+1; j < n; ++j){
                    coeff = pvector[ind];
                    A(i,j) = A(j,i) = coeff;
                    coeff = -vvector[ind];
                    B(i,j) = B(j,i) = coeff;
                    ++ind;
                }
            }
            precomputedValues.computed = true;
        }
    }

    NT positiveLinearIntersection(VT const & p, VT const & v) {
        createMatricesForPositiveLinearIntersection(p, v);
        return EigenvaluesProblem.minPosLinearEigenvalue(precomputedValues.B, precomputedValues.A,
                                                                precomputedValues.eigenvector);
    }

    /// Computes the intersection of the line a + tv,
    /// assuming we start at t=0 and that b has zero everywhere and 1 in its i-th coordinate.
    /// Solve the generalized eigenvalue problem A+tB = lmi(a + tv)
    /// So A = lmi(a) and B=lmi(v) - A0
    /// \param[in] a Input vector
    /// \param[in] coordinate Indicator of the i-th coordinate, 1 <= coordinate <= dimension
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT coordinateIntersection(VT const & a, int const coordinate) {
        // return EigenvaluesProblem.symGeneralizedProblem(matrix, *.getMatrix(coordinate)));
        return EigenvaluesProblem.symGeneralizedProblem(precomputedValues.A, *(lmi.getMatrix(coordinate)));
    }

    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int const& rand_coord,
                                          VT&)
    {
        return coordinateIntersection(r.getCoefficients(), rand_coord);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v)
    {
        NT pos_inter = positiveLinearIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: the eigenvector A(p)*e = 0
    /// \param[out] ret The unit normal vector at p
    void unit_normal(VT p, VT const& e, VT &ret) const
    {   
        NT sum_sqqrt_sq = NT(0);
        for (int ind = 0; ind < d; ++ind) {
            std::pair<int,int> indices = getMatrixIndices(n, ind);
            ret(i) = e[indices.first]*e[indices.second];
            sum_sqqrt_sq += ret(ind)*ret(ind);
        }
        ret /= std::sqrt(sum_sqqrt_sq); //normalize
    }

    /// Test if a point p is in the spectrahedron
    /// \param p is the current point
    /// \return true if position is outside the spectrahedron
    int is_in(Point const& p, NT tol=NT(0)) {  
        return !isExterior(p.getCoefficients());
    }

    bool isExterior(VT const & p) {
        if(precomputedValues.computed_A = false){
            buildMatrix(p.getCoefficients(), n, precomputedValues.A);
        }
        return isExterior(precomputedValues.A);
    }


    /*
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

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT&)
    {
        return line_positive_intersect(r, v);
    }

    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(PointType &r,
                                          PointType&,
                                          unsigned int const& rand_coord,
                                          unsigned int&,
                                          VT&)
    {
        return coordinateIntersection(r.getCoefficients(), rand_coord);
    }
    */

};