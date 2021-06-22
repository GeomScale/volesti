// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2021 Vaibhav Thakkar

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include "poset.h"
#include <Eigen/Eigen>
#include "preprocess/max_inscribed_ball.hpp"
#ifndef VOLESTIPY
    #include "lp_oracles/solve_lp.h"
#endif


template <typename Point>
class OrderPolytope {
public:
    typedef Point PointType;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    
private:
    Poset poset;
    unsigned int _d;    // dimension
    
    VT b;
    VT row_norms;
    MT _A;  // representing as Ax <= b for ComputeInnerBall and printing

    unsigned int _num_hyperplanes;
    bool _normalized;

public:
    OrderPolytope(const Poset& _poset) : poset(_poset)
    {
        _d = poset.num_elem();
        _num_hyperplanes = 2*_d + poset.num_relations(); // 2*d are for >=0 and <=1 constraints
        b = Eigen::MatrixXd::Zero(_num_hyperplanes, 1);
        row_norms = Eigen::MatrixXd::Constant(_num_hyperplanes, 1, 1.0);

        // first add (ai >= 0) or (-ai <= 0) rows    
        _A.topLeftCorner(_d, _d) = -Eigen::MatrixXd::Identity(_d, _d);

        // next add (ai <= 1) rows
        _A.block(_d, 0, _d, _d) = Eigen::MatrixXd::Identity(_d, _d);
        b.block(_d, 0, _d, 1) = Eigen::MatrixXd::Constant(_d, 1, 1.0);

        // next add the relations
        unsigned int num_relations = poset.num_relations();
        for(int idx=0; idx<num_relations; ++idx) {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(idx);
            _A(2*_d + idx, curr_relation.first)  = 1;
            _A(2*_d + idx, curr_relation.second) = -1;
        }
        row_norms.block(2*_d, 0, num_relations, 1) = Eigen::MatrixXd::Constant(num_relations, 1, sqrt(2));

        _normalized = false;
    }


    // return dimension
    unsigned int dimension() const
    {
        return _d;
    }


    // return number of hyperplanes
    unsigned int num_hyperplanes() const 
    {
        return _num_hyperplanes;
    }


    // print polytope in Ax <= b format
    void print()
    {
        std::cout << " " << _A.rows() << " " << _d << " double" << std::endl;
        for (unsigned int i = 0; i < _A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                if (!_normalized)
                    std::cout << _A(i, j) << " ";
                else
                    std::cout << _A(i, j) / row_norms(i) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
        }
    }


    /** multiply the sparse matrix A of the order polytope by a vector x
     *  if transpose = false     : return Ax
     *  else if transpose = true : return (A^T)x
     */
    VT vec_mult(const VT& x, bool transpose=false) {
        unsigned int rows = _num_of_hyperplanes;
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
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(i - 2*_d);

            if (!transpose) {
                if (! _normalized)
                    res(i) = x(curr_relation.first) - x(curr_relation.second);
                else
                    res(i) = (x(curr_relation.first) - x(curr_relation.second)) / row_norms(i);
            }
            else {
                if (! _normalized) {
                    res(curr_relation.first)  += x(i);
                    res(curr_relation.second) -= x(i);    
                }
                else {
                    res(curr_relation.first)  += x(i) / row_norms(i);
                    res(curr_relation.second) -= x(i) / row_norms(i);
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

        for (int i = 0; i < _d; i++) {
            // check violation of point between 0 and 1
            if (pt_coeffs(i) < NT(-tol) || (pt_coeffs(i) - 1.0) > NT(tol))
                return 0;
        }

        // check violations of order relations
        if (!poset.is_in(pt_coeffs, tol) {
            return 0;
        }

        return -1;
    }


    //Compute Chebyshev ball of the polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point, NT> ComputeInnerBall()
    {
       normalize();

        // change entries of A, doing here as won't be required in 
        // optimized volume calculation of order-polytope
        MT A = _A.rowwise().normalized();

        #ifndef VOLESTIPY   // as _A is never normalized in closed form
            _innerball = ComputeChebychevBall<NT, Point>(A, b); // use lpsolve library
        #else

            if (_innerball.second < 0.0) {

                NT const tol = 0.00000001;
                std::tuple<VT, NT, bool> innerball = max_inscribedball(A, b, 150, tol);

                // check if the solution is feasible
                if (is_in(Point(std::get<0>(innerball))) == 0 || std::get<1>(innerball) < NT(0) ||
                    std::isnan(std::get<1>(innerball)) || std::isinf(std::get<1>(innerball)) ||
                    !std::get<2>(innerball) || is_inner_point_nan_inf(std::get<0>(innerball)))
                {
                    _innerball.second = -1.0;
                } else
                {
                    _innerball.first = Point(std::get<0>(innerball));
                    _innerball.second = std::get<1>(innerball);
                }
            }
        #endif

        return _innerball;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the order polytope
    std::pair<NT,NT> line_intersect(Point const& r, Point const& v, bool pos = false) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom, sum_denom;
        
        unsigned int i = 0;
        int rows = _num_of_hyperplanes, facet;
        
        sum_nom.noalias() = b - vec_mult(r.getCoefficients());
        sum_denom.noalias() = vec_mult(v.getCoefficients());

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();

        // iterate over all hyperplanes
        for(; i<rows; ++i) {
            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
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


    // compute intersection point of ray starting from r and pointing to v
    // with the order-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v, const VT &Ar,
            const VT &Av) const 
    {
        return line_intersect(r, v);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the order-polytope
    std::pair<NT,NT> line_intersect(const Point &r, const Point &v, const VT &Ar,
            const VT &Av, const NT &lambda_prev) const         
    {
        return line_intersect(r, v);
    }


    std::pair<NT, int> line_positive_intersect(const Point &r, const Point &v, const VT &Ar,
                                               const VT &Av, const NT &lambda_prev) const {
        return line_positive_intersect(r, v);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev) const
    {
        return line_intersect(r, v, Ar, Av, lambda_prev, true);
    }


    //-------------------------accelarated billiard--------------------------------//
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters &params) const
    {
        return line_positive_intersect(r, v);
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
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters &params) const
    {
        return line_positive_intersect(r, v);
    }
    //------------------------------------------------------------------------------//


    // Compute the intersection of a coordinate ray
    // with the order polytope
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          VT& lamdas) const
    {
        Point v(_d);
        v.set_coord(rand_coord, 1.0);
        return line_intersect<NT>(r, v);
    }


    std::pair<NT,NT> line_intersect_coord(const Point &r,
                                          const Point &r_prev,
                                          const unsigned int rand_coord,
                                          const unsigned int rand_coord_prev,
                                          const VT &lamdas) const {
        return line_intersect_coord(r, rand_coord, lamdas);
    }



    // REMOVED: as this destroys sparsity
    // // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    // void linear_transformIt(MT const& T)
    // {
    //     _A = _A * T;
    // }

    
    // no points given for the rounding, you have to sample from the polytope
    template <typename T>
    bool get_points_for_rounding (T const& /*randPoints*/)
    {
        return false;
    }


    // shift polytope by a point c
    void shift(const VT &c)
    {
        b -= vec_mult(c);
    }


    void normalize()
    {
        if (_normalized = true)
            return;

        // for b and _A, first 2*_d rows are already normalized, for 
        _normalized = true; // -> will be used to make changes in entries of _A
        for (unsigned int i = 0; i < _num_hyperplanes; ++i)
        {
            b(i) = b(i) / row_norms(i);
        }
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(NT const& radius) const
    {        
        std::vector <NT> dists(_num_hyperplanes, NT(0));

        for (unsigned int i = 0; i < _num_hyperplanes; ++i) {
            if (_normalized)
                dists(i) = b(i);
            else
                dists(i) = b(i) / row_norms(i);
        }

        return dists;
    }


    // compute reflection in O(1) time for order polytope
    void compute_reflection(Point& v, Point const&, int const& facet) const
    {
        NT dot_prod;
        if (facet < _d) {
            dot_prod = -v[facet];
        }
        else if (facet < 2*_d)  {
            dot_prod = v[facet - _d];
        }
        else {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(facet - 2*_d);
            dot_prod = v[curr_relation.first] - v[curr_relation.second];
            dot_prod = dot_prod / row_norms(facet - 2*_d);
        }

        // calculating -> v += -2 * dot_prod * A.row(facet);
        if (facet < _d) {
            v.set_coord(facet, v[facet] - 2 * dot_prod * (-1.0));
        }
        else if (facet < 2*_d)  {
            v.set_coord(facet-_d, v[facet-_d] - 2 * dot_prod * (1.0));
        }
        else {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(facet - 2*_d);
            v.set_coord(curr_relation.first, v[curr_relation.first] - 2 * dot_prod * (1.0));
            v.set_coord(curr_relation.second, v[curr_relation.second] - 2 * dot_prod * (-1.0));
        }
    }


    template <typename update_parameters>
    void compute_reflection(Point &v, const Point &, update_parameters const& params) const 
    {
        NT dot_prod = params.inner_vi_ak;
        int facet = params.facet_prev;

        // calculating -> v += -2 * dot_prod * A.row(facet);
        if (facet < _d) {
            v.set_coord(facet, v[facet] - 2 * dot_prod * (-1.0));
        }
        else if (facet < 2*_d)  {
            v.set_coord(facet-_d, v[facet-_d] - 2 * dot_prod * (1.0));
        }
        else {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(facet - 2*_d);
            v.set_coord(curr_relation.first, v[curr_relation.first] - 2 * dot_prod * (1.0));
            v.set_coord(curr_relation.second, v[curr_relation.second] - 2 * dot_prod * (-1.0));
        }
    }


    NT log_barrier(Point &x, NT t = NT(100)) const
    {
      VT slack_vec = b - vec_mult(x.getCoefficients());
      NT total = slack_vec.log().sum();
                                        
      return total / t;
    }


    // calculated -ve grad log barrier
    Point grad_log_barrier(Point &x, NT t = NT(100)) 
    {
      VT slack_vec = b - vec_mult(x.getCoefficients());
      VT total_vec = vec_mult(slack_vec.inverse(), true);
      Point total(total_vec);

      total = (1.0 / t) * total;
      return total;
    }
};
