#ifndef VOLESTI_CORRE_SPECTRAHEDRON2_H
#define VOLESTI_CORRE_SPECTRAHEDRON2_H

#include "corre_matrix.hpp"
//
/// This class handles the spectrahedra of correlation matrices
/// @tparam Point
template<typename CorreMatrix>
class CorreSpectra2 : public Spectrahedron<CorreMatrix> {
    public:

    /// The numeric/matrix/vector types we use
    typedef CorreMatrix                                         PointType;
    typedef typename PointType::FT                              NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

    /// The size of the matrix
    unsigned int n;

    VT eigenvector;

    /// Constructor of correlation matrix spectrahedra
    /// \param[in] : matrix size
    CorreSpectra2(unsigned int n){
        int i,j;
        this->n = n;
        this->d = n*(n-1)/2;
        this->_inner_ball.first = PointType(this->d);
        this->_inner_ball.second = 1/std::sqrt(this->d);
        this->eigenvector.setZero(n);
    }

    /// \returns The size of the matrix
    unsigned int matrixSize() const {
        return n;
    }

    std::pair<PointType, NT> getInnerBall() const {
        return this->_inner_ball;
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] r A point on the boundary of the spectrahedron
    /// \param[in] v The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    template <typename update_parameters>
    void compute_reflection(PointType &v, PointType const& r, update_parameters& ) const {
        MT grad = MT::Zero(this->n, this->n);
        int i, j;
        NT sum_sq = NT(0), dot = NT(0);

        for(i = 0; i < n ; ++i){
            for(j = 0; j < i; ++j){
                grad(i,j) = eigenvector[i]*eigenvector[j];
                sum_sq += grad(i,j)*grad(i,j);
                dot += grad(i,j) * v.mat(i,j);
            }
        }
        dot = 2 * dot / sum_sq;
        
        // std::cout << "dot = " << dot << std::endl;
        // for(i = 0; i < n ; ++i){
        //     for(j = 0; j < i; ++j){
        //         NT tmp = dot*grad(i,j);
        //         // grad(i,j) =  tmp;
        //         std::cout << grad(i,j)  << std::endl;
        //     }
        // }
        
        grad = dot*grad;
        v -= PointType(grad);
    }

    NT positiveLinearIntersection(PointType const & r, PointType const & v) {
        return this->EigenvaluesProblem.minPosLinearEigenvalue2(r.mat, (-v).mat, eigenvector);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v) {
        NT pos_inter = positiveLinearIntersection(r, v);
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

    /// Test if a point p is in the spectrahedron
    /// \param p is the current point
    /// \return true if position is outside the spectrahedron
    int is_in(PointType const& p, NT tol=NT(0)) {  
        return !isExterior(p);
    }

    bool isExterior(MT const& mat) {
        return !this->EigenvaluesProblem.isPositiveSemidefinite(-mat);
    }
};

#endif //VOLESTI_CORRE_SPECTRAHEDRON_H