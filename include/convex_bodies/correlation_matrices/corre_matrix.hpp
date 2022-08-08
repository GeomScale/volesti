#ifndef VOLESTI_CORRE_MATRIX_H
#define VOLESTI_CORRE_MATRIX_H

/// This class handles the PointType used by CorreSpectra_MT class.
/// Every point is a correlation matrix and only the lower triangular part is stored.

/// @tparam NT
template<typename NT>
class CorreMatrix{
    public:

    /// The numeric/matrix/vector types we use
    typedef NT                                                  FT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic,1>                 Coeff;

    MT mat;

    CorreMatrix(){}

    CorreMatrix(unsigned int n){
        mat = MT::Identity(n,n);
    }
    
    CorreMatrix(MT const& mat){
        this->mat = mat;
    }

    void operator+= (const CorreMatrix<NT> & p)
    {
        this->mat += p.mat;
    }

    void operator-= (const CorreMatrix<NT> & p)
    {
        this->mat -= p.mat;
    }

    void operator= (const CorreMatrix<NT> & p)
    {
        this->mat = p.mat;
    }

    CorreMatrix<NT> operator- () const
    {
        CorreMatrix<NT> temp;
        temp.mat = - this->mat;
        return temp;
    }

    void operator*= (const FT k)
    {
        this->mat = k*this->mat;
    }

    CorreMatrix<NT> operator* (const FT k) const
    {
        MT M = k * this->mat;
        return CorreMatrix<NT>(M);
    }

    void operator/= (const FT k)
    {
        this->mat /= k;
    }

    NT dot(MT grad){
        int i, j, n = this->mat.rows();
        NT ret = NT(0);
        for(i = 0; i < n ; ++i){
            for(j = 0; j < i; ++j){
                ret += this->mat(i,j) * grad(i,j);
            }
        }
        return ret;
    }

    void print() const
    {
        int n = this->mat.rows(), i, j;
        for(i = 0; i < n ; ++i){
            for(j = 0; j < i; ++j){
                std::cout<< this->mat(i,j) <<" ";
            }
        }
        std::cout<<"\n";
    }

    VT getCoefficients(){
        int n = this->mat.rows(), ind = 0, i, j;
        VT coeff(n*(n-1)/2);
        for(i = 0; i < n ; ++i){
            for(j = 0; j < i; ++j){
                coeff(ind) = this->mat(i,j);
                ++ind;
            }
        }
        return coeff;
    }
};

template<typename NT>
CorreMatrix<NT> operator* (const NT k, CorreMatrix<NT> p)
{
    return p * k;
}

template<typename NT>
std::ostream& operator<<(std::ostream& os, const CorreMatrix<NT>& p)
{
    os << p.mat;
    return os;
}

#endif //VOLESTI_CORRE_MATRIX_H
