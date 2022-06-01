/// This class handles correlation matrices
/// @tparam NT Numeric Type

template<typename NT>
class Point_matrix {
    private:
    
        /// The size of the matrix
        unsigned int n;

        /// The length of coeffs
        unsigned int d;
    
    public:

    /// The numeric/matrix/vector types we use
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;

    /// The correlation matrix
    MT matrix;

    
    /// The coefficients of the matrix
    VT coeffs;

    /// Constructors

    Point_matrix(int n){
        this->n = n;
        d = n*(n-1)/2;
        coeffs = VT::Zero(d);
        matrix = MT::Identity(n,n);
    }

    Point_matrix(int n, VT coeffs){
        this->n = n;
        this->coeffs = coeffs;
        d = n*(n-1)/2;
        unsigned int i, j, ind;
        NT coeff;
        matrix = MT::Identity(n,n);
        for(i = 0; i < n ; ++i){
            for(j = i+1; j < n; ++j){
                ind = ((((n<<1)-i-2)*(i+1)) >> 1)  + j - n;
                coeff = coeffs(ind);
                matrix(i,j) = coeff;
                matrix(j,i) = coeff;
            }
        }
    }

    /// \returns The dimension of coefficient vector
    unsigned int dimension() const {
        return d;
    }

    /// \returns The size of the matrix
    unsigned int getSizeOfMatrix() const {
        return n;
    }

    /// \return The associated correlation matrix
    MT getMatrix() const {
        return matrix;
    }

    /// Prints the correlation matrix
    void print() const {
        std::cout << matrix << std::endl;
    }
};