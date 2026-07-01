#ifndef SIMPLEXINTERSECTBALL_H
#define SIMPLEXINTERSECTBALL_H

#include <limits>
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <Eigen/Eigen>

/// This class represents the intersection of a simplex with the unit ball.
/// The simplex is given in H-representation:
///     A x <= b
/// and the ball is:
///     ||x - x0|| <= 1
///
/// Return convention for is_in:
///     -1 : inside
///      0 : outside
///
/// \tparam Point Point type used by volesti
/// \tparam MT_type Matrix type for A
template
<
    typename Point,
    typename MT_type = Eigen::Matrix<typename Point::FT, Eigen::Dynamic, Eigen::Dynamic>
>
class SimplexIntersectBall {
public:
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef MT_type                                           MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;

        static NT pi()
    {
        return std::acos(NT(-1));
    }

private:
    unsigned int _d;  // dimension
    MT           A;   // matrix A
    VT           b;   // vector b, such that A x <= b
    MT           V;   // simplex vertices, stored column-wise
    VT           x0;  // center of the unit ball
    VT           Vnorms;

public:
    SimplexIntersectBall() {}

    SimplexIntersectBall(unsigned int d_,
                         MT const& A_,
                         VT const& b_,
                         MT const& V_,
                         VT const& x0_)
        : _d{d_}, A{A_}, b{b_}, V{V_}, x0{x0_}
    {
        Vnorms = V.colwise().norm();
        Vnorms = Vnorms.cwiseProduct(Vnorms);
    }

    // Copy constructor
    SimplexIntersectBall(SimplexIntersectBall<Point, MT_type> const& p)
        : _d{p._d}, A{p.A}, b{p.b}, V{p.V}, x0{p.x0}, Vnorms{p.Vnorms}
    {
    }

    // Return dimension
    unsigned int dimension() const
    {
        return _d;
    }

    // Return number of facets / hyperplanes
    int num_of_hyperplanes() const
    {
        return A.rows();
    }

    // Return matrix A
    MT get_mat() const
    {
        return A;
    }

    // Return vector b
    VT get_vec() const
    {
        return b;
    }

    // Return simplex vertices
    MT get_vertices() const
    {
        return V;
    }

    // Return center of the unit ball
    VT get_center() const
    {
        return x0;
    }

    // Change matrix A
    void set_mat(MT const& A2)
    {
        A = A2;
    }

    // Change vector b
    void set_vec(VT const& b2)
    {
        b = b2;
    }

    // Change simplex vertices
    void set_vertices(MT const& V2)
    {
        V = V2;
        Vnorms = V.colwise().norm();
        Vnorms = Vnorms.cwiseProduct(Vnorms);
    }

    // Change ball center
    void set_center(VT const& x02)
    {
        x0 = x02;
    }

    // Check if Point p lies in the simplex-ball intersection:
    //     A p <= b
    //     ||p - x0|| <= 1
    int is_in(Point const& p, NT tol = NT(0)) const
    {
        VT p_vec = p.getCoefficients();

        // Check ball condition
        VT diff = p_vec - x0;
        NT radius_tol = NT(1) + tol;
        if (diff.squaredNorm() > radius_tol * radius_tol) return 0;

        // Check simplex inequalities A p <= b
        VT temp = b - A * p_vec;
        const NT* Ax_b_data = temp.data();

        for (int i = 0; i < A.rows(); i++) {
            if ((*Ax_b_data) < NT(-tol)) {
                return 0;
            }
            Ax_b_data++;
        }

        return -1;
    }

    int is_in_optimized(Point const& p,
                         VT& Ar,
                         VT& Av,
                         NT const& lambda_prev,
                         NT tol = NT(0)) const
    {
        VT p_vec = p.getCoefficients();

        // Ball membership check.
        VT diff = p_vec - x0;
        NT radius_tol = NT(1) + tol;
        if (diff.squaredNorm() > radius_tol * radius_tol) return 0;

        // Update cached A*x along the great-circle rotation:
        // x(lambda) = cos(lambda) * r + sin(lambda) * v.
        //
        // Here Ar stores A*r and Av stores A*v from the previous step.
        Ar.noalias() = std::cos(lambda_prev) * Ar + std::sin(lambda_prev) * Av;

        VT temp = b - Ar;
        const NT* Ax_b_data = temp.data();

        for (int i = 0; i < A.rows(); i++) {
            if ((*Ax_b_data) < NT(-tol)) return 0;
            Ax_b_data++;
        }

        return -1;
    }


    // Compute intersection parameters of the line r + lambda * v
    // with the simplex A x <= b.
    //
    // Returns:
    //   first  = smallest positive lambda
    //   second = largest negative lambda
    std::pair<NT, NT> line_intersect(Point const& r, Point const& v) const
    {
        NT lambda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        VT sum_nom;
        VT sum_denom;

        int m = num_of_hyperplanes();

        sum_nom.noalias() = b - A * r.getCoefficients();
        sum_denom.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {
            if (*sum_denom_data != NT(0)) {
                lambda = *sum_nom_data / *sum_denom_data;

                if (lambda < min_plus && lambda > NT(0)) {
                    min_plus = lambda;
                }

                if (lambda > max_minus && lambda < NT(0)) {
                    max_minus = lambda;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        return std::make_pair(min_plus, max_minus);
    }

    // Optimized version of line_intersect.
    // It also computes and returns:
    //     Ar = A * r
    //     Av = A * v
    //
    // If pos = false:
    //   returns (min positive lambda, max negative lambda)
    //
    // If pos = true:
    //   returns (min positive lambda, facet index)
    std::pair<NT, NT> line_intersect(Point const& r,
                                     Point const& v,
                                     VT& Ar,
                                     VT& Av,
                                     bool pos = false) const
    {
        NT lambda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        VT sum_nom;
        int m = num_of_hyperplanes();
        int facet = -1;

        Ar.noalias() = A * r.getCoefficients();
        sum_nom.noalias() = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data != NT(0)) {
                lambda = *sum_nom_data / *Av_data;

                if (lambda < min_plus && lambda > NT(0)) {
                    min_plus = lambda;
                    if (pos) {
                        facet = i;
                    }
                } else if (lambda > max_minus && lambda < NT(0)) {
                    max_minus = lambda;
                }
            }

            Av_data++;
            sum_nom_data++;
        }

        if (pos) {
            return std::make_pair(min_plus, facet);
        }

        return std::make_pair(min_plus, max_minus);
    }

    // Optimized line_intersect version using the previous lambda.
    // Ar is updated as:
    //     Ar <- Ar + lambda_prev * Av
    //
    // Then Av is recomputed as:
    //     Av <- A * v
    std::pair<NT, NT> line_intersect(Point const& r,
                                     Point const& v,
                                     VT& Ar,
                                     VT& Av,
                                     NT const& lambda_prev,
                                     bool pos = false) const
    {   
        (void)r;  
        NT lambda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        VT sum_nom;
        int m = num_of_hyperplanes();
        int facet = -1;

        Ar.noalias() += lambda_prev * Av;
        sum_nom.noalias() = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data != NT(0)) {
                lambda = *sum_nom_data / *Av_data;

                if (lambda < min_plus && lambda > NT(0)) {
                    min_plus = lambda;
                    if (pos) {
                        facet = i;
                    }
                } else if (lambda > max_minus && lambda < NT(0)) {
                    max_minus = lambda;
                }
            }

            Av_data++;
            sum_nom_data++;
        }

        if (pos) {
            return std::make_pair(min_plus, facet);
        }

        return std::make_pair(min_plus, max_minus);
    }

    // Compute intersection angles of the great circle
    // x(lambda) = cos(lambda) * r + sin(lambda) * v
    // with the simplex boundary A x <= b.
    //
    // Input:
    //   Ar = A * r
    //   Av = A * v
    //
    // Returns:
    //   first  = smallest positive angle lambda
    //   second = largest negative angle lambda
    std::pair<NT, NT> compute_intersections(VT& Ar, VT& Av) const
    {
        NT D;
        NT C1;
        NT C2;
        NT eval;

        NT max_root = std::numeric_limits<NT>::lowest();
        NT min_root = std::numeric_limits<NT>::max();

        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        int m = num_of_hyperplanes();

        bool set_negative_root = false;
        bool set_positive_root = false;
        bool pos_D = false;

        NT* Av_data = Av.data();
        NT* Ar_data = Ar.data();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            D = (*Ar_data) * (*Ar_data)
              + (*Av_data) * (*Av_data)
              - (*b_data) * (*b_data);

            if (D > NT(0)) {
                pos_D = true;

                NT denom = (*Ar_data) * (*Ar_data)
                         + (*Av_data) * (*Av_data);

                C1 = std::asin((((*Av_data) * (*b_data)) + ((*Ar_data) * std::sqrt(D))) / denom);
                C2 = std::asin((((*Av_data) * (*b_data)) - ((*Ar_data) * std::sqrt(D))) / denom);

                eval = (*Ar_data) * std::cos(C1) + (*Av_data) * std::sin(C1) - (*b_data);
                if (!(eval > -NT(1e-05) && eval < NT(1e-05))) {
                    C1 = pi() - C1;
                }

                if (C1 < min_plus && C1 > NT(0)) {
                    min_plus = C1;
                    set_positive_root = true;
                } else if (C1 > max_minus && C1 < NT(0)) {
                    max_minus = C1;
                    set_negative_root = true;
                }

                if (C1 > max_root && C1 < NT(2) * pi()) {
                    max_root = C1;
                }

                if ((C1 < min_root) && (C1 > (-NT(2) * pi()))) {
                    min_root = C1;
                }

                eval = (*Ar_data) * std::cos(C2) + (*Av_data) * std::sin(C2) - (*b_data);
                if (!(eval > -NT(1e-05) && eval < NT(1e-05))) {
                    C2 = pi() - C2;
                }

                if (C2 < min_plus && C2 > NT(0)) {
                    min_plus = C2;
                    set_positive_root = true;
                } else if (C2 > max_minus && C2 < NT(0)) {
                    max_minus = C2;
                    set_negative_root = true;
                }

                if (C2 > max_root && C2 < NT(2) * pi()) {
                    max_root = C2;
                }

                if (C2 < min_root && C2 > (-NT(2) * pi())) {
                    min_root = C2;
                }
            }

            Av_data++;
            Ar_data++;
            b_data++;
        }

        if (!set_negative_root) {
            if (pos_D) {
                max_minus = max_root - NT(2) * pi();
            } else {
                max_minus = NT(0);
            }
        }

        if (!set_positive_root) {
            if (pos_D) {
                min_plus = min_root + NT(2) * pi();
            } else {
                min_plus = NT(2) * pi();
            }
        }

        return std::make_pair(min_plus, max_minus);
    }

    // Great-circle intersection.
    // Computes Ar = A*r and Av = A*v, then calls compute_intersections.
    std::pair<NT, NT> gc_intersect(Point const& r,
                                   Point const& v,
                                   VT& Ar,
                                   VT& Av) const
    {
        Ar.noalias() = A * r.getCoefficients();
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections(Ar, Av);
    }

    // Optimized great-circle intersection using previous lambda.
    std::pair<NT, NT> gc_intersect(Point const& r,
                                   Point const& v,
                                   VT& Ar,
                                   VT& Av,
                                   NT const& lambda_prev) const
    {
        (void)r;

        Ar.noalias() = std::cos(lambda_prev) * Ar + std::sin(lambda_prev) * Av;
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections(Ar, Av);
    }

    // Great-circle intersection assuming Ar is already up to date.
    std::pair<NT, NT> gc_intersect_optimized(Point const& r,
                                             Point const& v,
                                             VT& Ar,
                                             VT& Av,
                                             NT const& lambda_prev) const
    {
        (void)r;
        (void)lambda_prev;

        Av.noalias() = A * v.getCoefficients();

        return compute_intersections(Ar, Av);
    }

    // Compute the first positive intersection angle of the great circle
    // x(lambda) = cos(lambda) * r + sin(lambda) * v
    // with the simplex boundary A x <= b.
    //
    // Input:
    //   Ar = A * r
    //   Av = A * v
    //
    // Returns:
    //   first  = smallest positive angle lambda
    //   second = index of the facet hit
    std::pair<NT, int> compute_intersections_positive(VT const& Ar, VT const& Av) const
    {
        NT D;
        NT C1;
        NT C2;
        NT eval;
        NT min_root = std::numeric_limits<NT>::max();
        NT min_plus = std::numeric_limits<NT>::max();

        int m = num_of_hyperplanes();
        int facet = -1;
        int facet_min = -1;

        bool set_positive_root = false;
        bool pos_D = false;

        const NT* Av_data = Av.data();
        const NT* Ar_data = Ar.data();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            D = (*Ar_data) * (*Ar_data)
              + (*Av_data) * (*Av_data)
              - (*b_data) * (*b_data);

            if (D > NT(0)) {
                pos_D = true;

                NT denom = (*Ar_data) * (*Ar_data)
                         + (*Av_data) * (*Av_data);

                C1 = std::asin((((*Av_data) * (*b_data)) + ((*Ar_data) * std::sqrt(D))) / denom);
                C2 = std::asin((((*Av_data) * (*b_data)) - ((*Ar_data) * std::sqrt(D))) / denom);

                eval = (*Ar_data) * std::cos(C1) + (*Av_data) * std::sin(C1) - (*b_data);
                if (!(eval > -NT(1e-05) && eval < NT(1e-05))) {
                    C1 = pi() - C1;
                }

                if (C1 < min_plus && C1 > NT(0)) {
                    min_plus = C1;
                    set_positive_root = true;
                    facet = i;
                }

                if ((C1 < min_root) && (C1 > (-NT(2) * pi()))) {
                    min_root = C1;
                    facet_min = i;
                }

                eval = (*Ar_data) * std::cos(C2) + (*Av_data) * std::sin(C2) - (*b_data);
                if (!(eval > -NT(1e-05) && eval < NT(1e-05))) {
                    C2 = pi() - C2;
                }

                if (C2 < min_plus && C2 > NT(0)) {
                    min_plus = C2;
                    set_positive_root = true;
                    facet = i;
                }

                if (C2 < min_root && C2 > (-NT(2) * pi())) {
                    min_root = C2;
                    facet_min = i;
                }
            }

            Av_data++;
            Ar_data++;
            b_data++;
        }

        if (!set_positive_root) {
            if (pos_D) {
                min_plus = min_root + NT(2) * pi();
                facet = facet_min;
            } else {
                min_plus = NT(2) * pi();
            }
        }

        return std::make_pair(min_plus, facet);
    }

    // Great-circle first positive intersection.
    // Computes Ar = A*r and Av = A*v, then calls compute_intersections_positive.
    std::pair<NT, int> gc_intersect_positive(Point const& r,
                                             Point const& v,
                                             VT& Ar,
                                             VT& Av) const
    {
        Ar.noalias() = A * r.getCoefficients();
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections_positive(Ar, Av);
    }

    // Optimized great-circle first positive intersection using previous lambda.
    std::pair<NT, int> gc_intersect_positive(Point const& r,
                                             Point const& v,
                                             VT& Ar,
                                             VT& Av,
                                             NT const& lambda_prev) const
    {
        (void)r;

        Ar.noalias() = std::cos(lambda_prev) * Ar + std::sin(lambda_prev) * Av;
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections_positive(Ar, Av);
    }


        // Compute all intersection roots of the great circle
    // x(lambda) = cos(lambda) * r + sin(lambda) * v
    // with the simplex boundary A x <= b.
    //
    // Returns:
    //   first  = negative roots in (-pi, 0)
    //   second = positive roots in (0, pi)
    std::pair<VT, VT> compute_intersections_all_roots(VT& Ar, VT& Av) const
    {
        std::vector<NT> neg_roots;
        std::vector<NT> pos_roots;

        int m = num_of_hyperplanes();

        NT* Av_data = Av.data();
        NT* Ar_data = Ar.data();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            NT denom = (*Ar_data) * (*Ar_data) + (*Av_data) * (*Av_data);
            NT D = denom - (*b_data) * (*b_data);

            if (D > NT(0)) {
                NT sqrtD = std::sqrt(D);

                NT C1 = std::asin((((*Av_data) * (*b_data)) + ((*Ar_data) * sqrtD)) / denom);
                NT C2 = std::asin((((*Av_data) * (*b_data)) - ((*Ar_data) * sqrtD)) / denom);

                NT eval = (*Ar_data) * std::cos(C1)
                        + (*Av_data) * std::sin(C1)
                        - (*b_data);

                if (!(eval > -NT(1e-05) && eval < NT(1e-05))) {
                    C1 = pi() - C1;
                }

                if (C1 > pi()) {
                    C1 -= NT(2) * pi();
                } else if (C1 < -pi()) {
                    C1 += NT(2) * pi();
                }

                if ((C1 > -pi()) && (C1 < NT(0))) {
                    neg_roots.push_back(C1);
                } else if ((C1 < pi()) && (C1 > NT(0))) {
                    pos_roots.push_back(C1);
                }

                eval = (*Ar_data) * std::cos(C2)
                     + (*Av_data) * std::sin(C2)
                     - (*b_data);

                if (!(eval > -NT(1e-05) && eval < NT(1e-05))) {
                    C2 = pi() - C2;
                }

                if (C2 > pi()) {
                    C2 -= NT(2) * pi();
                } else if (C2 < -pi()) {
                    C2 += NT(2) * pi();
                }

                if ((C2 > -pi()) && (C2 < NT(0))) {
                    neg_roots.push_back(C2);
                } else if ((C2 < pi()) && (C2 > NT(0))) {
                    pos_roots.push_back(C2);
                }
            }

            Av_data++;
            Ar_data++;
            b_data++;
        }

        std::sort(neg_roots.begin(), neg_roots.end());
        std::sort(pos_roots.begin(), pos_roots.end());

        VT neg_roots_vec(neg_roots.size());
        for (unsigned int i = 0; i < neg_roots.size(); i++) {
            neg_roots_vec(i) = neg_roots[i];
        }

        VT pos_roots_vec(pos_roots.size());
        for (unsigned int i = 0; i < pos_roots.size(); i++) {
            pos_roots_vec(i) = pos_roots[i];
        }

        return std::make_pair(neg_roots_vec, pos_roots_vec);
    }

    // Great-circle all-roots intersection.
    std::pair<VT, VT> gc_intersect_all_roots(Point const& r,
                                             Point const& v,
                                             VT& Ar,
                                             VT& Av) const
    {
        Ar.noalias() = A * r.getCoefficients();
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections_all_roots(Ar, Av);
    }

    // Optimized great-circle all-roots intersection using previous lambda.
    std::pair<VT, VT> gc_intersect_all_roots(Point const& r,
                                             Point const& v,
                                             VT& Ar,
                                             VT& Av,
                                             NT const& lambda_prev) const
    {
        (void)r;

        Ar.noalias() = std::cos(lambda_prev) * Ar + std::sin(lambda_prev) * Av;
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections_all_roots(Ar, Av);
    }

    // Apply linear transformation T to the simplex inequalities.
    // If the point transformation is x <- T^{-1} x,
    // then the H-representation matrix changes as A <- A * T.
    void linear_transformIt(MT const& T)
    {
        A = A * T;
    }

    // Shift the simplex by a vector c.
    // For A x <= b, after shifting x <- x + c,
    // the right-hand side changes as b <- b - A*c.
    void shift(VT const& c)
    {
        b -= A * c;
    }

    // Reflect tangent direction v at point p on the sphere
    // against the facet indexed by facet.
    //
    // p_p is the tangent-space projector:
    //     p_p = I - p p^T
    void compute_reflection(VT& v,
                            VT const& p,
                            MT const& p_p,
                            int const& facet) const
    {
        (void)p;

        VT u = p_p * A.row(facet).transpose();
        u *= (NT(1) / u.norm());

        v += -NT(2) * v.dot(u) * u;
    }

};
#endif