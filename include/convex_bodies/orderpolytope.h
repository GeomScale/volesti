#include <iostream>
#include "poset.h"


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
    
    VT _b;
    MT _A;  // representing as Ax <= b for ComputeInnerBall and printing

    unsigned int _num_hyperplanes;

public:
    OrderPolytope(const Poset& _poset) : poset(_poset)
    {
        _d = poset.num_elem();
        _num_hyperplanes = 2*_d + poset.num_relations(); // 2*d are for >=0 and <=1 constraints
        _b = Eigen::MatrixXd::Zero(_num_hyperplanes, 1);

        // first add (ai >= 0) or (-ai <= 0) rows    
        _A.topLeftCorner(_d, _d) = -Eigen::MatrixXd::Identity(_d, _d);

        // next add (ai <= 1) rows
        _A.block(_d, 0, _d, _d) = Eigen::MatrixXd::Identity(_d, _d);
        _b.block(_d, 0, _d, _d) = Eigen::MatrixXd::Constant(_d, 1, 1.0);

        // next add the relations
        for(int idx=0; idx<poset.num_relations(); ++idx) {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(idx);
            _A(2*_d + idx, curr_relation.first)  = 1;
            _A(2*_d + idx, curr_relation.second) = -1;
        }
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
                std::cout << _A(i, j) << " ";
            }
            std::cout << "<= " << _b(i) << std::endl;
        }
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


    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point, NT> ComputeInnerBall()
    {
       normalize();
        #ifndef VOLESTIPY
            _inner_ball = ComputeChebychevBall<NT, Point>(_A, _b); // use lpsolve library
        #else

            if (_inner_ball.second < 0.0) {

                NT const tol = 0.00000001;
                std::tuple<VT, NT, bool> inner_ball = max_inscribed_ball(_A, _b, 150, tol);

                // check if the solution is feasible
                if (is_in(Point(std::get<0>(inner_ball))) == 0 || std::get<1>(inner_ball) < NT(0) ||
                    std::isnan(std::get<1>(inner_ball)) || std::isinf(std::get<1>(inner_ball)) ||
                    !std::get<2>(inner_ball) || is_inner_point_nan_inf(std::get<0>(inner_ball)))
                {
                    _inner_ball.second = -1.0;
                } else
                {
                    _inner_ball.first = Point(std::get<0>(inner_ball));
                    _inner_ball.second = std::get<1>(inner_ball);
                }
            }
        #endif

        return _inner_ball;
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
        VT r_coeffs = r.getCoefficients();
        VT v_coeffs = v.getCoefficients();
        NT sum_nom_data, sum_denom_data;

        // first _d rows of >=0 constraints
        for(; i < _d; ++i) {
            sum_nom_data   = b(i) - (-1.0)*r_coeffs(i);
            sum_denom_data = (-1.0)*v_coeffs(i);

            if (sum_denom_data != NT(0)) {
                lamda = sum_nom_data / sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }
        }

        // next _d rows of <=1 constraints
        for(; i < 2*_d; ++i) {
            sum_nom_data   = b(i) - (1.0)*r_coeffs(i - _d);
            sum_denom_data = (1.0)*v_coeffs(i - _d);

            if (sum_denom_data != NT(0)) {
                lamda = sum_nom_data / sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }
        }

        // next rows are for order relations
        for(; i < rows; ++i) {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(i - 2*_d);
            
            sum_nom_data = b(i) - ((1.0)*r_coeffs(curr_relation.first) - (1.0)*r_coeffs(curr_relation.second));
            sum_denom_data = (1.0)*v_coeffs(curr_relation.first) - (1.0)*v_coeffs(curr_relation.second);

            if (sum_denom_data != NT(0)) {
                lamda = sum_nom_data / sum_denom_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }
                else if (lamda > max_minus && lamda < 0) {
                    max_minus = lamda;
                }
            }
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

};
