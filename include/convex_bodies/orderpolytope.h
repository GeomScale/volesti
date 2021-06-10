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
    MT _A;
    unsigned int _num_hyperplanes;

public:
    OrderPolytope(const Poset& _poset) : poset(_poset)
    {
        _d = poset.num_elem();
        _num_hyperplanes = 2*_d + poset.num_relations(); // 2*d are for >=0 and <=1 constraints
        _b = Eigen::MatrixXd::Zero(_num_hyperplanes, 1);
        _b.block(n, 0, n, 1) = Eigen::MatrixXd::Constant(n, 1, 1.0);
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
        unsigned int rows = _num_hyperplanes;
        std::cout << " " << rows << " " << _d << " double" << std::endl;
        unsigned int i = 0;

        // first _d rows of >=0 constraints
        for(; i < _d; ++i) {
            for (unsigned int j = 0; j < _d; ++j) {
                if (j == i)
                    std::cout << -1 << " ";
                else
                    std::cout << 0 << " ";
            }
            std::cout << "<= " << _b(i) << std::endl;
        }

        // next _d rows of <=1 constraints
        for(; i < 2*_d; ++i) {
            for (unsigned int j = 0; j < _d; ++j) {
                if ((j+_d) == i)
                    std::cout << 1 << " ";
                else
                    std::cout << 0 << " ";
            }
            std::cout << "<= " << _b(i) << std::endl;
        }

        // next rows are for order relations
        for(; i < rows; ++i) {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(i - 2*_d);
            for (unsigned int j = 0; j < _d; ++j) {
                if (j == curr_relation.first)
                    std::cout << 1 << " ";
                else if (j == curr_relation.second)
                    std::cout << -1 << " ";
                else
                    std::cout << 0 << " ";
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


    // compute intersection point of ray starting from r and pointing to v
    // with the order polytope
    std::pair<NT,NT> line_intersect(Point const& r, Point const& v) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom, sum_denom;
        
        unsigned int i = 0;
        int rows = _num_of_hyperplanes;
        VT r_coeffs = r.getCoefficients();
        VT v_coeffs = v.getCoefficients();
        NT sum_nom_data, sum_denom_data;

        // first _d rows of >=0 constraints
        for(; i < _d; ++i) {
            sum_nom_data   = b(i) - (-1.0)*r_coeffs(i);
            sum_denom_data = (-1.0)*v_coeffs(i);

            if (sum_denom_data != NT(0)) {
                lamda = sum_nom_data / sum_denom_data;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }

        // next _d rows of <=1 constraints
        for(; i < 2*_d; ++i) {
            sum_nom_data   = b(i) - (1.0)*r_coeffs(i - _d);
            sum_denom_data = (1.0)*v_coeffs(i - _d);

            if (sum_denom_data != NT(0)) {
                lamda = sum_nom_data / sum_denom_data;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }

        // next rows are for order relations
        for(; i < rows; ++i) {
            std::pair<unsigned int, unsigned int> curr_relation = poset.get_relation(i - 2*_d);
            
            sum_nom_data = b(i) - ((1.0)*r_coeffs(curr_relation.first) - (1.0)*r_coeffs(curr_relation.second));
            sum_denom_data = (1.0)*v_coeffs(curr_relation.first) - (1.0)*v_coeffs(curr_relation.second);

            if (sum_denom_data != NT(0)) {
                lamda = sum_nom_data / sum_denom_data;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }

        return std::make_pair(min_plus, max_minus);
    }

};
