/// This class handles points corresponding to correlation matrices
/// @tparam Point Type
template<typename K>
class Point_matrix : public Point<K>
{
    private:

        /// The numeric/matrix/vector types we use

        typedef typename K::FT 	FT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;

        /// The coefficients of the matrix
        Eigen::Matrix<, Eigen::Dynamic,1> coeffs;
        /// The correlation matrix
        MT matrix;

        /// The gradient vector
        VT grad;

        typedef typename std::vector<typename K::FT>::iterator iter;

        /// The size of the matrix
        unsigned int n;

        /// The length of coeffs
        unsigned int d;
    
    public:

    /// Constructors

    Point_matrix() {}

    Point_matrix(const unsigned int dim) {
        d = dim;
        coeffs.setZero(d);
    }

    point(const unsigned int dim, iter begin, iter endit)
    {
        d = dim;
        coeffs.resize(d);
        int i = 0;

        for (iter it=begin ; it != endit ; it++)
            coeffs(i++) = *it;
    }

    point(const Coeff& coeffs)
    {
            d = coeffs.rows();
            this->coeffs = coeffs;
    }

    point(const unsigned int dim, std::vector<typename K::FT> cofs)
    {
        d = dim;
        coeffs.resize(d);
        iter it = cofs.begin();
        int i=0;

        for (; it != cofs.end(); it++,i++)
            coeffs(i) = *it;

    }

    Point_matrix(int n){
        this->n = n;
        d = n*(n-1)/2;
        coeffs = Point(d);
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

void add(const Coeff& coeffs)
    {
        this->coeffs += coeffs;
    }

    const Coeff& getCoefficients() const
    {
        return coeffs;
    }

    /// \returns The dimension of coefficient vector
    int dimension() const {
        return d;
    }

    void set_dimension(const unsigned int dim) {
        d = dim;
    }

    void set_coord(const unsigned int i, FT coord) {
        coeffs(i) = coord;
    }

    void set_coeffs (const Coeff& coeffs2) {
        d = coeffs2.rows();
        coeffs = coeffs2;
    }

    void set_to_origin() {
        coeffs.setZero(d);
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

    
    
    
    FT operator[] (const unsigned int i) const
    {
        return coeffs(i);
    }

    FT* pointerToData()
    {
        return coeffs.data();
    }

    FT sum() const {
        return coeffs.sum();
    }

    void operator+= (const point& p)
    {
        coeffs += p.getCoefficients();
    }

    void operator+= (const Coeff& coeffs)
    {
        this->coeffs = coeffs + this->coeffs;
    }

    void operator= (const Coeff& coeffs)
    {
        this->coeffs = coeffs;
        d = coeffs.rows();
    }

    //TODO: avoid point construction in operators +,-,*
    point operator+ (const point& p) const
    {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs + p.getCoefficients();
        return temp;
    }

    point operator- (const point& p) const
    {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs - p.getCoefficients();
        return temp;
    }

    point operator* (const FT k) const
    {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs * k;
        return temp;
    }

    void operator*= (const FT k)
    {
        coeffs *= k;
    }

    void operator/= (const FT k)
    {
        coeffs /= k;
    }

    bool operator== (point& p) const
    {
        int i=0;
        const Coeff & coeffs = p.getCoefficients();

        /* degree of approximation in
        "The art of computer programming" (vol II), p. 234, Donald. E. Knuth. */
        FT e = 0.00000000001;
        for (i=0 ; i<d ; i++) {
            if (std::abs(this->coeffs(i) - coeffs(i)) > e *std::abs(this->coeffs(i)) ||
                    std::abs(this->coeffs(i) - coeffs(i)) > e *std::abs(coeffs(i)))
                return false;
        }

        return true;
    }

    FT distance(point const & p) {
        return (this->coeffs - p.coeffs).norm();
    }

    FT dot(const point& p) const
    {
        return coeffs.dot(p.getCoefficients());
    }

    FT dot(const Coeff& coeffs) const
    {
        return this->coeffs.dot(coeffs);
    }

    FT squared_length() const {
        FT lsq = length();
        return lsq * lsq;
    }

    FT length() const {
        return coeffs.norm();
    }

    void print() const
    {
        for(unsigned int i=0; i<d; i++){
            std::cout<<coeffs(i)<<" ";
        }
        std::cout<<"\n";
    }

    static point all_ones(int dim) {
      point p(dim);
      for (int i = 0; i < dim; i++) p.set_coord(i, 1.0);
      return p;
    }

};

template<typename K>
point<K> operator* (const typename K::FT& k, point<K> const& p)
{
    return p * k;
}

#endif

};