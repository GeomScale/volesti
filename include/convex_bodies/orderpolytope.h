// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2021 Vaibhav Thakkar

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ORDER_POLYTOPE_H
#define ORDER_POLYTOPE_H

#include <iostream>
#include "misc/poset.h"
#include <Eigen/Eigen>
#include "preprocess/max_inscribed_ball.hpp"
#ifndef DISABLE_LPSOLVE
    #include "lp_oracles/solve_lp.h"
#endif

/// This class represents an order polytope parameterized by a point type
/// \tparam Point Point type
template <typename Point>
class OrderPolytope {
public:
    typedef Point PointType;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

private:
    Poset _poset;
    unsigned int _d;    // dimension

    VT b;
    VT _row_norms;
    MT _A;  // representing as Ax <= b for ComputeInnerBall and printing

    unsigned int _num_hyperplanes;
    bool _normalized;

public:
    OrderPolytope(Poset const& poset) : _poset(poset)
    {
        _d = _poset.num_elem();
        _num_hyperplanes = 2*_d + _poset.num_relations(); // 2*d are for >=0 and <=1 constraints
        b = Eigen::MatrixXd::Zero(_num_hyperplanes, 1);
        _A = Eigen::MatrixXd::Zero(_num_hyperplanes, _d);
        _row_norms = Eigen::MatrixXd::Constant(_num_hyperplanes, 1, 1.0);

        // first add (ai >= 0) or (-ai <= 0) rows
        _A.topLeftCorner(_d, _d) = -Eigen::MatrixXd::Identity(_d, _d);

        // next add (ai <= 1) rows
        _A.block(_d, 0, _d, _d) = Eigen::MatrixXd::Identity(_d, _d);
        b.block(_d, 0, _d, 1) = Eigen::MatrixXd::Constant(_d, 1, 1.0);

        // next add the relations
        unsigned int num_relations = _poset.num_relations();
        for(int idx=0; idx<num_relations; ++idx) {
            std::pair<unsigned int, unsigned int> curr_relation = _poset.get_relation(idx);
            _A(2*_d + idx, curr_relation.first)  = 1;
            _A(2*_d + idx, curr_relation.second) = -1;
        }
        _row_norms.block(2*_d, 0, num_relations, 1) = Eigen::MatrixXd::Constant(num_relations, 1, sqrt(2));

        _normalized = false;
    }


    // return dimension
    unsigned int dimension() const
    {
        return _d;
    }


    // return number of hyperplanes
    unsigned int num_of_hyperplanes() const
    {
        return _num_hyperplanes;
    }


    // get ith column of A
    VT get_col (unsigned int i) const {
        return _A.col(i);
    }


    Eigen::SparseMatrix<NT> get_mat() const
    {
        return _A.sparseView();
    }


    VT get_vec() const
    {
        return b;
    }


    // print polytope in Ax <= b format
    void print() const
    {
        std::cout << " " << _A.rows() << " " << _d << " double" << std::endl;
        for (unsigned int i = 0; i < _A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                std::cout << _A(i, j) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
        }
    }


    /** multiply the sparse matrix A of the order polytope by a vector x
     *  if transpose = false     : return Ax
     *  else if transpose = true : return (A^T)x
     */
    VT vec_mult(VT const& x, bool transpose=false) const
    {
        unsigned int rows = num_of_hyperplanes();
        unsigned int i = 0;
        VT res;
        if (!transpose) res = Eigen::MatrixXd::Zero(rows, 1);
        else            res = Eigen::MatrixXd::Zero(_d, 1);

        // ------- no effect of normalize on first 2*_d rows, norm = 1 ---------
        // first _d rows of >=0 constraints
        for(; i < _d; ++i) {
            res(i) = (-1.0) * x(i);
        }

        // next _d rows of <=1 constraints
        for(; i < 2*_d; ++i) {
            if (!transpose) {
                res(i) = (1.0) * x(i - _d);
            }
            else {
                res(i - _d) += (1.0) * x(i);
            }
        }
        // -----------------------------------------------------------------

        // next rows are for order relations
        for(; i < rows; ++i) {
            std::pair<unsigned int, unsigned int> curr_relation = _poset.get_relation(i - 2*_d);

            if (!transpose) {
                if (! _normalized)
                    res(i) = x(curr_relation.first) - x(curr_relation.second);
                else
                    res(i) = (x(curr_relation.first) - x(curr_relation.second)) / _row_norms(i);
            }
            else {
                if (! _normalized) {
                    res(curr_relation.first)  += x(i);
                    res(curr_relation.second) -= x(i);
                }
                else {
                    res(curr_relation.first)  += x(i) / _row_norms(i);
                    res(curr_relation.second) -= x(i) / _row_norms(i);
                }
            }
        }

        return res;
    }


    //Check if Point p is in the order-polytope
    int is_in(Point const& p, NT tol=NT(0)) const
    {
        assert(p.dimension() == _d);

        VT pt_coeffs = p.getCoefficients();
        NT diff;

        for (int i = 0; i < _d; i++) {
            // DON'T JUST check violation of point between 0 and 1
            // as b will change for shifted polytope.
            diff = -pt_coeffs(i) - b(i);
            if (diff > NT(tol)) return 0;

            diff = pt_coeffs(i) - b(i + _d);
            if (diff > NT(tol)) return 0;
        }

        // check violations of order relations
        // again note that b can be arbitrary because of shifting
        unsigned int num_relations = _poset.num_relations();
        for(int idx=0; idx<num_relations; ++idx) {
            std::pair<unsigned int, unsigned int> curr_relation = _poset.get_relation(idx);
            diff = (pt_coeffs(curr_relation.first) - pt_coeffs(curr_relation.second));

            if (_normalized) diff /= _row_norms(idx + 2*_d);
            if((diff - b(idx + 2*_d)) > NT(tol))
                return 0;
        }

        return -1;
    }


    //Compute Chebyshev ball of the polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point, NT> ComputeInnerBall()
    {
        normalize();
        std::pair<Point, NT> inner_ball;
        #ifndef DISABLE_LPSOLVE
            inner_ball = ComputeChebychevBall<NT, Point>(_A, b); // use lpsolve library
        #else

            if (inner_ball.second <= NT(0)) {

                NT const tol = 0.00000001;
                std::tuple<VT, NT, bool> inner_ball = max_inscribed_ball(_A, b, 150, tol);

                // check if the solution is feasible
                if (is_in(Point(std::get<0>(inner_ball))) == 0 || std::get<1>(inner_ball) < NT(0) ||
                    std::isnan(std::get<1>(inner_ball)) || std::isinf(std::get<1>(inner_ball)) ||
                    !std::get<2>(inner_ball) || is_inner_point_nan_inf(std::get<0>(inner_ball)))
                {
                    inner_ball.second = -1.0;
                } else
                {
                    inner_ball.first = Point(std::get<0>(inner_ball));
                    inner_ball.second = std::get<1>(inner_ball);
                }
            }
        #endif

        return inner_ball;

    }


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    // compute intersection point of ray starting from r and pointing to v
    // with the order polytope
    std::pair<NT,NT> line_intersect(Point const& r, Point const& v, bool pos = false) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom, sum_denom;

        int rows = num_of_hyperplanes(), facet;

        sum_nom.noalias() = b - vec_mult(r.getCoefficients());
        sum_denom.noalias() = vec_mult(v.getCoefficients());

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        if (pos)
            return std::make_pair(min_plus, facet);

        return std::make_pair(min_plus, max_minus);
    }


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    // compute intersection point of ray starting from r and pointing to v
    // with the order-polytope
    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT& Ar,
                                    VT& Av,
                                    bool pos = false) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom;

        int rows = num_of_hyperplanes(), facet;

        Ar.noalias() = vec_mult(r.getCoefficients());
        Av.noalias() = vec_mult(v.getCoefficients());

        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = Av.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        if (pos)
            return std::make_pair(min_plus, facet);

        return std::make_pair(min_plus, max_minus);
    }


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    // compute intersection point of ray starting from r and pointing to v
    // with the order-polytope
    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT &Ar,
                                    VT &Av,
                                    NT const& lambda_prev,
                                    bool pos = false) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom;

        int rows = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        Av.noalias() = vec_mult(v.getCoefficients());

        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = Av.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        if (pos)
            return std::make_pair(min_plus, facet);

        return std::make_pair(min_plus, max_minus);
    }


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    // compute intersection point of a ray starting from r and pointing to v
    // with the order-polytope
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av) const
    {
        return line_intersect(r, v, Ar, Av, true);
    }


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    // compute intersection point of a ray starting from r and pointing to v
    // with the order-polytope
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev) const
    {
        return line_intersect(r, v, Ar, Av, lambda_prev, true);
    }


    //-------------------------accelerated billiard--------------------------------//
    // compute intersection point of a ray starting from r and pointing to v
    // with the order-polytope
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters &params) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        VT sum_nom;

        int rows = num_of_hyperplanes(), facet;

        Ar.noalias() = vec_mult(r.getCoefficients());
        Av.noalias() = vec_mult(v.getCoefficients());

        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = Av.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    params.inner_vi_ak = *sum_denom_data;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        params.facet_prev = facet;
        return std::make_pair(min_plus, facet);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               MT const& AA,
                                               update_parameters &params) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        VT sum_nom;

        int rows = num_of_hyperplanes(), facet;
        NT inner_prev = params.inner_vi_ak;

        Ar.noalias() += lambda_prev*Av;
        if(params.hit_ball) {
            Av.noalias() += (-2.0 * inner_prev) * (Ar / params.ball_inner_norm);
        } else {
            Av.noalias() += (-2.0 * inner_prev) * AA.col(params.facet_prev);
        }
        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = Av.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    params.inner_vi_ak = *sum_denom_data;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        params.facet_prev = facet;
        return std::make_pair(min_plus, facet);
    }


    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters &params) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        VT sum_nom;

        int rows = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        Av.noalias() = vec_mult(v.getCoefficients());

        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = Av.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    params.inner_vi_ak = *sum_denom_data;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        params.facet_prev = facet;
        return std::make_pair(min_plus, facet);
    }
    //------------------------------------------------------------------------------//


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    // Compute the intersection of a coordinate ray
    // with the order polytope
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          VT& lamdas) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_denom;

        int rows = num_of_hyperplanes();

        sum_denom = get_col(rand_coord);
        lamdas = b - vec_mult(r.getCoefficients());

        NT* sum_nom_data = lamdas.data();
        NT* sum_denom_data = sum_denom.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
        }

        return std::make_pair(min_plus, max_minus);
    }


    // TODO: This can be removed as only modified Accelerated Billiard Walk will be used with Order Polytope
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          VT& lamdas) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        int rows = num_of_hyperplanes();

        lamdas.noalias() += get_col(rand_coord_prev)
                        * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
        NT* sum_nom_data = lamdas.data();

        // iterate over all hyperplanes
        for(unsigned int i = 0; i<rows; ++i) {
            NT a = _A(i, rand_coord);
            if(_normalized) {
                a = a / _row_norms(i);
            }

            if (a == NT(0)) {
                // throw std::runtime_error("Error: division by 0");
            } else {
                lamda = *sum_nom_data / a;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }

            sum_nom_data++;
        }

        return std::make_pair(min_plus, max_minus);
    }


    // no points given for the rounding, you have to sample from the polytope
    template <typename T>
    bool get_points_for_rounding (T const& /*randPoints*/)
    {
        return false;
    }


    // shift polytope by a point c
    void shift(VT const& c)
    {
        b -= vec_mult(c);
    }


    // get a point inside the order polytope,
    // NOTE: the current implementation only works for non shifted order polytope
    VT inner_point()
    {
        // get topologically sorted list of indices
        std::vector<unsigned int> sorted_list = _poset.topologically_sorted_list();

        // vector to hold n linearly spaced values between 0-1
        std::vector<NT> lin_space_values(_d);
        NT start = 0.05, end = 0.95;
        NT h = (end - start)/static_cast<NT>(_d-1);
        NT val = start;
        for(auto x=lin_space_values.begin(); x!=lin_space_values.end(); ++x) {
            *x = val;
            val += h;
        }

        // final result vector
        VT res(_d);
        for(int i=0; i<_d; ++i) {
            unsigned int curr_idx = sorted_list[i];
            res(curr_idx) = lin_space_values[i];
        }

        return res;
    }


    void normalize()
    {
        if (_normalized == true)
            return;

        // for b and _A, first 2*_d rows are already normalized, for
        _normalized = true; // -> will be used to make normalization idempotent
        for (unsigned int i = 0; i < _num_hyperplanes; ++i)
        {
            _A.row(i) /= _row_norms(i);
            b(i) /= _row_norms(i);
        }
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(NT const& radius) const
    {
        std::vector <NT> dists(_num_hyperplanes, NT(0));

        for (unsigned int i = 0; i < _num_hyperplanes; ++i) {
            if (_normalized)
                dists[i] = b(i);
            else
                dists[i] = b(i) / _row_norms(i);
        }

        return dists;
    }


    // compute reflection given dot product and facet
    void compute_reflection(Point& v, NT dot_prod, unsigned int facet) const
    {
        // calculating -> v += -2 * dot_prod * A.row(facet);
        if (facet < _d) {
            v.set_coord(facet, v[facet] - 2 * dot_prod * (-1.0));
        }
        else if (facet < 2*_d)  {
            v.set_coord(facet-_d, v[facet-_d] - 2 * dot_prod * (1.0));
        }
        else {
            std::pair<unsigned int, unsigned int> curr_relation = _poset.get_relation(facet - 2*_d);
            v.set_coord(curr_relation.first, v[curr_relation.first] - 2 * dot_prod * (1.0 / _row_norms(facet)));
            v.set_coord(curr_relation.second, v[curr_relation.second] - 2 * dot_prod * (-1.0 / _row_norms(facet)));
        }
    }


    // compute reflection in O(1) time for order polytope
    void compute_reflection(Point& v, Point const&, unsigned int facet) const
    {
        NT dot_prod;
        if (facet < _d) {
            dot_prod = -v[facet];
        }
        else if (facet < 2*_d)  {
            dot_prod = v[facet - _d];
        }
        else {
            std::pair<unsigned int, unsigned int> curr_relation = _poset.get_relation(facet - 2*_d);
            dot_prod = v[curr_relation.first] - v[curr_relation.second];
            dot_prod = dot_prod / _row_norms(facet);
        }

        compute_reflection(v, dot_prod, facet);
    }


    template <typename update_parameters>
    void compute_reflection(Point &v, Point const&, update_parameters const& params) const
    {
        NT dot_prod = params.inner_vi_ak;
        unsigned int facet = params.facet_prev;
        compute_reflection(v, dot_prod, facet);
    }
};

#endif
