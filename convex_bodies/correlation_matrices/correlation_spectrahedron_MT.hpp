// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_CONVEX_BODIES_CORRELATION_MATRICES_VOLESTI_CORRELATION_SPECTRAHEDRON_MT_HPP
#define VOLESTI_CONVEX_BODIES_CORRELATION_MATRICES_VOLESTI_CORRELATION_SPECTRAHEDRON_MT_HPP

#include "convex_bodies/correlation_matrices/corre_matrix.hpp"

/// This class handles the spectrahedra of correlation matrices
/// @tparam CorreMatrix The Correlation Matrix
template<typename CorreMatrix>
class CorrelationSpectrahedron_MT : public Spectrahedron<CorreMatrix>{
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

    CorrelationSpectrahedron_MT(unsigned int n){
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
    void compute_reflection(PointType &v, PointType const &r, update_parameters&) const {
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
        grad = dot*grad;
        v -= PointType(grad);
    }

    /// Computes the minimal positive t s.t. r+t*v intersects the boundary of the spectrahedron
    /// \param[in] r
    /// \param[in] v
    /// \param[out] a NT value t
    NT positiveLinearIntersection(PointType const &r, PointType const &v){

        // minPosLinearEigenvalue_EigenSymSolver(A,B) computes the minimal positive eigenvalue of A-t*B

        return this->EigenvaluesProblem.minPosLinearEigenvalue_EigenSymSolver(r.mat, (-v).mat, eigenvector);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const &r,
                                               PointType const &v) {
        NT pos_inter = positiveLinearIntersection(r, v);
        return std::pair<NT, int> (pos_inter, -1);
    }

    std::pair<NT, int> line_positive_intersect(PointType const &r,
                                               PointType const &v,
                                               VT&,
                                               VT& ,
                                               NT const&){
        return line_positive_intersect(r, v);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const &r,
                                               PointType const &v,
                                               VT&,
                                               VT&){
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const &r,
                                               PointType const &v,
                                               VT&,
                                               VT& ,
                                               NT const&,
                                               update_parameters&){
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const &r,
                                               PointType const &v,
                                               VT&,
                                               VT&,
                                               NT const&,
                                               MT const&,
                                               update_parameters&){
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType const &r,
                                                     PointType const &v,
                                                     VT&,
                                                     VT&,
                                                     update_parameters&){
        return line_positive_intersect(r, v);
    }

    // compute intersection point of ray starting from r and pointing to v
    std::pair<NT,NT> line_intersect(PointType const &r, PointType const &v) const {
        return this->EigenvaluesProblem.symGeneralizedProblem(-r.mat, -v.mat);
    }

    std::pair<NT,NT> line_intersect(PointType const &r,
                                    PointType const &v,
                                    VT&,
                                    VT&) const {
        return line_intersect(r, v);
    }

    std::pair<NT,NT> line_intersect(PointType const &r,
                                    PointType const &v,
                                    VT&,
                                    VT&,
                                    NT&) const {
        return line_intersect(r, v);
    }


    /// Test if a point p is in the spectrahedron
    /// \param p is the current point
    /// \return true if position is outside the spectrahedron
    int is_in(PointType const &p, NT tol=NT(0)) const {
        if(this->EigenvaluesProblem.isPositiveSemidefinite(p.mat)){
            return -1;
        }
        return 0;
    }

    bool isExterior(MT const &mat) const {
        return !this->EigenvaluesProblem.isPositiveSemidefinite(mat);
    }

    MT get_mat() const {
        return MT::Identity(this->d, this->d);
    }
};

#endif //VOLESTI_CONVEX_BODIES_CORRELATION_MATRICES_VOLESTI_CORRELATION_SPECTRAHEDRON_MT_HPP
