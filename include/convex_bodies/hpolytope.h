// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018-19 programs.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.
//Contributed and/or modified by Luca Perju, as part of Google Summer of Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HPOLYTOPE_H
#define HPOLYTOPE_H

#include <limits>
#include <iostream>
#include <Eigen/Eigen>
#include "preprocess/max_inscribed_ball.hpp"
#include "root_finders/quadratic_polynomial_solvers.hpp"
#ifndef DISABLE_LPSOLVE
    #include "lp_oracles/solve_lp.h"
#endif



// check if an Eigen vector contains NaN or infinite values
template <typename VT>
bool is_inner_point_nan_inf(VT const& p)
{
    typedef Eigen::Array<bool, Eigen::Dynamic, 1> VTint;
    VTint a = p.array().isNaN();
    for (int i = 0; i < p.rows(); i++) {
        if (a(i) || std::isinf(p(i))){
            return true;
        }
    }
    return false;
}

/// This class describes a polytope in H-representation or an H-polytope
/// i.e. a polytope defined by a set of linear inequalities
/// \tparam Point Point type
template 
<
    typename Point, 
    typename MT_type = Eigen::Matrix<typename Point::FT, Eigen::Dynamic, Eigen::Dynamic>
>
class HPolytope {
public:
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef typename std::vector<NT>::iterator                viterator;
    typedef MT_type                                           MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;

private:
    unsigned int         _d; //dimension
    MT                   A; //matrix A
    VT                   b; // vector b, s.t.: Ax<=b
    std::pair<Point, NT> _inner_ball;
    bool                 normalized = false; // true if the polytope is normalized
    bool                 has_ball = false;

public:
    //TODO: the default implementation of the Big3 should be ok. Recheck.
    HPolytope() {}

    HPolytope(unsigned d_, MT const& A_, VT const& b_) :
        _d{d_}, A{A_}, b{b_}
    {
    }

    template<typename T = DenseMT>
    HPolytope(unsigned d_, DenseMT const& A_, VT const& b_, typename std::enable_if<!std::is_same<MT, T>::value, T>::type* = 0) :
        _d{d_}, A{A_.sparseView()}, b{b_}
    {
    }

    // Copy constructor
    HPolytope(HPolytope<Point, MT> const& p) :
            _d{p._d}, A{p.A}, b{p.b}, _inner_ball{p._inner_ball}, normalized{p.normalized}, has_ball{p.has_ball}
    {
    }

    //define matrix A and vector b, s.t. Ax<=b,
    // from a matrix that contains both A and b, i.e., [A | b ]
    HPolytope(std::vector<std::vector<NT>> const& Pin)
    {
        _d = Pin[0][1] - 1;
        A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                A.coeffRef(i - 1, j - 1) = -Pin[i][j];
            }
        }
        has_ball = false;
        //_inner_ball = ComputeChebychevBall<NT, Point>(A, b);
    }


    std::pair<Point, NT> InnerBall() const
    {
        return _inner_ball;
    }

    void set_InnerBall(std::pair<Point,NT> const& innerball) //const
    {
        _inner_ball = innerball;
        has_ball = true;
    }

    void set_interior_point(Point const& r)
    {
        _inner_ball.first = r;
        if(!is_in(r)) {
            std::cerr << "point outside of polytope was set as interior point" << std::endl;
            has_ball = false;
            return;
        }
        if(is_normalized()) {
            _inner_ball.second = (b - A * r.getCoefficients()).minCoeff();
        } else {
            _inner_ball.second = std::numeric_limits<NT>::max();
            for(int i = 0; i < num_of_hyperplanes(); ++i) {
                NT dist = (b(i) - A.row(i).dot(r.getCoefficients()) ) / A.row(i).norm();
                if(dist < _inner_ball.second) {
                    _inner_ball.second = dist;
                }
            }
        }
        has_ball = true;
    }

    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //First try using max_inscribed_ball
    //Use LpSolve library if it fails
    std::pair<Point, NT> ComputeInnerBall()
    {
        normalize();
        if (!has_ball) {
            
            has_ball = true;
            NT const tol = 1e-08;
            std::tuple<VT, NT, bool> inner_ball = max_inscribed_ball(A, b, 5000, tol);

            // check if the solution is feasible
            if (is_in(Point(std::get<0>(inner_ball))) == 0 || std::get<1>(inner_ball) < tol/2.0 ||
                std::isnan(std::get<1>(inner_ball)) || std::isinf(std::get<1>(inner_ball)) ||
                is_inner_point_nan_inf(std::get<0>(inner_ball))) {
                
                std::cerr << "Failed to compute max inscribed ball, trying to use lpsolve" << std::endl;
                #ifndef DISABLE_LPSOLVE
                    _inner_ball = ComputeChebychevBall<NT, Point>(A, b); // use lpsolve library
                #else
                    std::cerr << "lpsolve is disabled, unable to compute inner ball";
                    has_ball = false;
                #endif
            } else {
                _inner_ball.first = Point(std::get<0>(inner_ball));
                _inner_ball.second = std::get<1>(inner_ball);
            }
        }
        return _inner_ball;
    }

    // return dimension
    unsigned int dimension() const
    {
        return _d;
    }


    // return the number of facets
    int num_of_hyperplanes() const
    {
        return A.rows();
    }

    int num_of_generators() const
    {
        return 0;
    }


    // return the matrix A
    MT get_mat() const
    {
        return A;
    }


    MT get_AA() const {
        return A * A.transpose();
    }

    // return the vector b
    VT get_vec() const
    {
        return b;
    }

    bool is_normalized ()
    {
        return normalized;
    }


    // change the matrix A
    void set_mat(MT const& A2)
    {
        A = A2;
        normalized = false;
        has_ball = false;
    }


    // change the vector b
    void set_vec(VT const& b2)
    {
        b = b2;
        has_ball = false;
    }

    Point get_mean_of_vertices() const
    {
        return Point(_d);
    }

    NT get_max_vert_norm() const
    {
        return 0.0;
    }


    // print polytope in input format
    void print() {
        std::cout << " " << A.rows() << " " << _d << " double" << std::endl;
        for (unsigned int i = 0; i < A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                std::cout << A.coeff(i, j) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
        }
    }


    // Compute the reduced row echelon form
    // used to transofm {Ax=b,x>=0} to {A'x'<=b'}
    // e.g. Birkhoff polytopes
    /*
    // Todo: change the implementation in order to use eigen matrix and vector.
    int rref(){
        to_reduced_row_echelon_form(_A);
        std::vector<int> zeros(_d+1,0);
        std::vector<int> ones(_d+1,0);
        std::vector<int> zerorow(_A.size(),0);
        for (int i = 0; i < _A.size(); ++i)
        {
            for (int j = 0; j < _d+1; ++j){
                if ( _A[i][j] == double(0)){
                    ++zeros[j];
                    ++zerorow[i];
                }
                if ( _A[i][j] == double(1)){
                    ++ones[j];
                }
            }
        }
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            int j =0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ){
                if(zeros[j]==_A.size()-1 && ones[j]==1)
                    (*mit).erase(lit);
                else{ //reverse sign in all but the first column
                    if(lit!=mit->end()-1) *lit = (-1)*(*lit);
                    ++lit;
                }
                ++j;
            }
        }
        //swap last and first columns
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            double temp=*(mit->begin());
            *(mit->begin())=*(mit->end()-1);
            *(mit->end()-1)=temp;
        }
        //delete zero rows
        for (typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ){
            int zero=0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit){
                if(*lit==double(0)) ++zero;
            }
            if(zero==(*mit).size())
                _A.erase(mit);
            else
                ++mit;
        }
        //update _d
        _d=(_A[0]).size();
        // add unit vectors
        for(int i=1;i<_d;++i){
            std::vector<double> e(_d,0);
            e[i]=1;
            _A.push_back(e);
        }
        // _d should equals the dimension
        _d=_d-1;
        return 1;
    }*/


    //Check if Point p is in H-polytope P:= Ax<=b
    int is_in(Point const& p, NT tol=NT(0)) const
    {
        int m = A.rows();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            //Check if corresponding hyperplane is violated
            if (*b_data - A.row(i) * p.getCoefficients() < NT(-tol))
                return 0;

            b_data++;
        }
        return -1;
    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point const& r, Point const& v) const
    {

        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom, sum_denom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes();


        sum_nom.noalias() = b - A * r.getCoefficients();
        sum_denom.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {

            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }

            sum_nom_data++;
            sum_denom_data++;
        }
        return std::make_pair(min_plus, max_minus);
    }

    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
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
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() = A * r.getCoefficients();
        sum_nom = b - Ar;
        Av.noalias() = A * v.getCoefficients();;


        NT* Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }

            Av_data++;
            sum_nom_data++;
        }
        if (pos) return std::make_pair(min_plus, facet);
        return std::make_pair(min_plus, max_minus);
    }

    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT& Ar,
                                    VT& Av,
                                    NT const& lambda_prev,
                                    bool pos = false) const
    {

        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT  sum_nom;
        NT mult;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        sum_nom = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
            Av_data++;
            sum_nom_data++;
        }
        if (pos) return std::make_pair(min_plus, facet);
        return std::make_pair(min_plus, max_minus);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av) const
    {
        return line_intersect(r, v, Ar, Av, true);
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


    //---------------------------accelarated billiard----------------------------------
    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters& params) const
    {
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        NT lamda = 0;
        VT sum_nom;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() = A * r.getCoefficients();
        sum_nom.noalias() = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    params.inner_vi_ak = *Av_data;
                }
            }

            Av_data++;
            sum_nom_data++;
        }
        params.facet_prev = facet;
        return std::pair<NT, int>(min_plus, facet);
    }


    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     NT const& lambda_prev,
                                                     DenseMT const& AA,
                                                     update_parameters& params) const
    {

        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        NT lamda = 0;
        NT inner_prev = params.inner_vi_ak;
        VT sum_nom;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        if(params.hit_ball) {
            Av.noalias() += (-2.0 * inner_prev) * (Ar / params.ball_inner_norm);
        } else {
            Av.noalias() += ((-2.0 * inner_prev) * AA.col(params.facet_prev));
        }
        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    params.inner_vi_ak = *Av_data;
                }
            }
            Av_data++;
            sum_nom_data++;
        }
        params.facet_prev = facet;
        return std::pair<NT, int>(min_plus, facet);
    }



    template <typename update_parameters, typename set_type, typename AA_type>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                                     VT& Ar,
                                                     VT& Av,
                                                     NT const& lambda_prev,
                                                     set_type& distances_set,
                                                     AA_type const& AA,
                                                     update_parameters& params) const
    {
        NT inner_prev = params.inner_vi_ak;
        
        NT* Av_data = Av.data();
        const NT* b_data = b.data();

        for (Eigen::SparseMatrix<double>::InnerIterator it(AA, params.facet_prev); it; ++it) {

            *(Av_data + it.row()) += (-2.0 * inner_prev) * it.value();
            NT val = distances_set.get_val(it.row());
            val += (val - params.moved_dist) * 2.0 * inner_prev * it.value() / *(Av_data + it.row());
            distances_set.change_val(it.row(), val, params.moved_dist);
        }
        
        std::pair<NT, int> ans = distances_set.get_min();
        ans.first -= params.moved_dist;

        params.inner_vi_ak = *(Av_data + ans.second);
        params.facet_prev = ans.second;

        /*if(ans.first < 0.00000001) {
            std::cout << "distance of 0 found" << std::endl;
            exit(0);
        }*/

        return ans;
    }


    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters& params) const
    {
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        NT lamda = 0;
        VT sum_nom;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        sum_nom.noalias() = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    params.inner_vi_ak = *Av_data;
                }
            }
            Av_data++;
            sum_nom_data++;
        }
        params.facet_prev = facet;
        return std::pair<NT, int>(min_plus, facet);
    }

    //-----------------------------------------------------------------------------------//


    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          VT& lamdas) const
    {

        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_denom;

        int m = num_of_hyperplanes();

        sum_denom = A.col(rand_coord);
        lamdas.noalias() = b - A * r.getCoefficients();

        NT* lamda_data = lamdas.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {

            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<sum_denom<<std::endl;
                ;
            } else {
                lamda = *lamda_data * (1 / *sum_denom_data);
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            lamda_data++;
            sum_denom_data++;
        }
        return std::make_pair(min_plus, max_minus);
    }


    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          VT& lamdas) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();

        int m = num_of_hyperplanes();

        lamdas.noalias() += (DenseMT)(A.col(rand_coord_prev)
                         * (r_prev[rand_coord_prev] - r[rand_coord_prev]));
        NT* data = lamdas.data();

        for (int i = 0; i < m; i++) {
            NT a = A.coeff(i, rand_coord);

            if (a == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *data / a;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            data++;
        }
        return std::make_pair(min_plus, max_minus);
    }


    //------------------------------oracles for exponential sampling---------------//////

    std::pair<NT, int> get_positive_quadratic_root(Point const& r, //current poistion
                                                   Point const& v, // current velocity
                                                   VT& Ac, // the product Ac where c is the bias vector of the exponential distribution
                                                   NT const& T, // the variance of the exponential distribution
                                                   VT& Ar, // the product Ar
                                                   VT& Av, // the product Av
                                                   int& facet_prev) const //the facet that the trajectory hit in the previous reflection
    {
        NT lamda = 0;
        NT lamda2 =0;
        NT lamda1 =0;
        NT alpha;
        NT min_plus  = std::numeric_limits<NT>::max();
        VT sum_nom;
        int m = num_of_hyperplanes();
        int facet = -1;

        sum_nom = Ar - b;
        Av.noalias() = A * v.getCoefficients();;

        NT* Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();
        NT* Ac_data = Ac.data();

        for (int i = 0; i < m; i++)
        {
            alpha = -((*Ac_data) / (2.0 * T));
            if (solve_quadratic_polynomial(alpha, (*Av_data), (*sum_nom_data), lamda1, lamda2))
            {
                lamda = pick_first_intersection_time_with_boundary(lamda1, lamda2, i, facet_prev);
                if (lamda < min_plus && lamda > 0)
                {
                    min_plus = lamda;
                    facet = i;
                }
            }
            Av_data++;
            sum_nom_data++;
            Ac_data++;
        }
        facet_prev = facet;
        return std::make_pair(min_plus, facet);
    }


    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> quadratic_positive_intersect(Point const& r, //current poistion
                                                    Point const& v, // current velocity
                                                    VT& Ac, // the product Ac where c is the bias vector of the exponential distribution
                                                    NT const& T, // the variance of the exponential distribution
                                                    VT& Ar, // the product Ar
                                                    VT& Av, // the product Av
                                                    int& facet_prev) const //the facet that the trajectory hit in the previous reflection
    {
        Ar.noalias() = A * r.getCoefficients();
        return get_positive_quadratic_root(r, v, Ac, T, Ar, Av, facet_prev);
    }

    std::pair<NT, int> quadratic_positive_intersect(Point const& r, //current poistion
                                                    Point const& v, // current velocity
                                                    VT& Ac,  // the product Ac where c is the bias vector of the exponential distribution
                                                    NT const& T, // the variance of the exponential distribution
                                                    VT& Ar, // the product Ar
                                                    VT& Av, // the product Av
                                                    NT const& lambda_prev, // the intersection time of the previous reflection
                                                    int& facet_prev) const //the facet that the trajectory hit in the previous reflection
    {
        Ar.noalias() += ((lambda_prev * lambda_prev) / (-2.0*T)) * Ac + lambda_prev * Av;
        return get_positive_quadratic_root(r, v, Ac, T, Ar, Av, facet_prev);
    }

    NT pick_first_intersection_time_with_boundary(NT const& lamda1, NT const& lamda2, int const& current_facet, int const& previous_facet) const
    {
        if (lamda1 == lamda2)
        {
            return lamda1;
        }
        NT lamda;
        const double tol = 1e-10;
        std::pair<NT, NT> minmax_values = std::minmax(lamda1, lamda2);

        lamda = (previous_facet == current_facet)
            ? minmax_values.second < NT(tol) ? minmax_values.first : minmax_values.second
            : minmax_values.second;

        if (lamda1 * lamda2 < NT(0))
        {
            lamda = (previous_facet == current_facet)
            ? (minmax_values.second < NT(tol)) ? minmax_values.first : minmax_values.second
            : minmax_values.second;
        }
        else
        {
            lamda = (previous_facet == current_facet)
            ? (minmax_values.first >= NT(0) && minmax_values.first < NT(tol))
            ? minmax_values.second : minmax_values.first
            : minmax_values.first;
        }
        return lamda;
    }


    // Boundary oracle for exact hmc spherical gaussian sampling
    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b (the ray describes a curve)
    std::pair<NT, int> trigonometric_positive_intersect(Point const& r, Point const& v,
                                                        NT const& omega, int &facet_prev) const
    {
        constexpr NT pi_2 = NT(2.0) * M_PI;
        NT t = std::numeric_limits<NT>::max();

        int m = num_of_hyperplanes();
        int facet = -1;

        VT sum_nom;
        VT sum_denom;
        sum_nom.noalias() = A * r.getCoefficients();
        sum_denom.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();
        const NT* b_data = b.data();

        const NT omega_sqr = omega * omega;
        const NT pi_2_omega = pi_2 / omega;

        for (int i = 0; i < m; i++) {

            NT C = std::sqrt((*sum_nom_data) * (*sum_nom_data) + ((*sum_denom_data) * (*sum_denom_data))
              / omega_sqr);
            NT Phi = std::atan((-(*sum_denom_data)) / ((*sum_nom_data) * omega));

            if ((*sum_denom_data) < 0.0 && Phi < 0.0) {
                Phi += M_PI;
            } else if ((*sum_denom_data) > 0.0 && Phi > 0.0) {
                Phi -= M_PI;
            }

            if (C > (*b_data)) {
                NT acos_b = std::acos((*b_data) / C);
                NT t1 = (acos_b - Phi) / omega;
                if (facet_prev == i && std::abs(t1) < 1e-10){
                    t1 = pi_2_omega;
                }

                NT t2 = (-acos_b - Phi) / omega;
                if (facet_prev == i && std::abs(t2) < 1e-10){
                    t2 = pi_2_omega;
                }

                t1 += (t1 < NT(0)) ? pi_2_omega : NT(0);
                t2 += (t2 < NT(0)) ? pi_2_omega : NT(0);

                NT tmin = std::min(t1, t2);

                if (tmin < t && tmin > NT(0)) {
                    facet = i;
                    t = tmin;
                }
            }

            sum_nom_data++;
            sum_denom_data++;
            b_data++;
        }
        facet_prev = facet;
        return std::make_pair(t, facet);
    }


    // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    template<typename T_type>
    void linear_transformIt(T_type const& T)
    {
        if constexpr (std::is_same<MT, DenseMT>::value) {
            A = A * T;
        } else {
            A = (A * T).sparseView();
        }
        normalized = false;
        has_ball = false;
    }


    // shift polytope by a point c

    void shift(const VT &c)
    {
        b -= A*c;
        has_ball = false;
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(NT const& radius) const
    {
        unsigned int i=0;
        std::vector <NT> dists(num_of_hyperplanes(), NT(0));
        typename std::vector<NT>::iterator disit = dists.begin();
        for ( ; disit!=dists.end(); disit++, i++)
            *disit = b(i) / A.row(i).norm();

        return dists;
    }

    // no points given for the rounding, you have to sample from the polytope
    template <typename T>
    bool get_points_for_rounding (T const& /*randPoints*/)
    {
        return false;
    }

    MT get_T() const
    {
        return A;
    }

    void normalize()
    {
        if(normalized)
            return;
        NT row_norm;
        for (int i = 0; i < A.rows(); ++i) {
            row_norm = A.row(i).norm();
            if (row_norm != 0.0) {
                A.row(i) /= row_norm;
                b(i) /= row_norm;
            }
        }
        normalized = true;
    }

    void compute_reflection(Point& v, Point const&, int const& facet) const
    {
        v += -2 * v.dot(A.row(facet)) * A.row(facet);
    }

    void resetFlags() {}

    NT log_barrier(Point &x, NT t = NT(100)) const {
      int m = num_of_hyperplanes();
      NT total = NT(0);
      NT slack;

      for (int i = 0; i < m; i++) {
        slack = b(i) - x.dot(A.row(i));
        total += log(slack);
      }

      return total / t;
    }

    Point grad_log_barrier(Point &x, NT t = NT(100)) {
      int m = num_of_hyperplanes();
      NT slack;

      Point total(x.dimension());

      for (int i = 0; i < m; i++) {
        slack = b(i) - x.dot(A.row(i));
        total = total + (1 / slack) * A.row(i);
      }
      total = (1.0 / t) * total;
      return total;
    }

    template <typename update_parameters>
    void compute_reflection(Point &v, Point const&, update_parameters const& params) const {
            Point a((-2.0 * params.inner_vi_ak) * A.row(params.facet_prev));
            v += a;
    }

    template <typename update_parameters>
    auto compute_reflection(Point &v, Point &p, update_parameters const& params) const
         -> std::enable_if_t<std::is_same_v<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>> && !std::is_same_v<update_parameters, int>, void> { // MT must be in RowMajor format
            NT* v_data = v.pointerToData();
            NT* p_data = p.pointerToData();
            for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, params.facet_prev); it; ++it) {
                *(v_data + it.col()) += (-2.0 * params.inner_vi_ak) * it.value();
                *(p_data + it.col()) -= (-2.0 * params.inner_vi_ak * params.moved_dist) * it.value();
            }
    }

    template <typename update_parameters>
    NT compute_reflection(Point &v, const Point &, DenseMT const &AE, VT const &AEA, NT const &vEv, update_parameters const &params) const {

            Point a((-2.0 * params.inner_vi_ak) * A.row(params.facet_prev));
            VT x = v.getCoefficients();
            NT new_vEv = vEv - (4.0 * params.inner_vi_ak) * (AE.row(params.facet_prev).dot(x) - params.inner_vi_ak * AEA(params.facet_prev));
            v += a;
            NT coeff = std::sqrt(vEv / new_vEv);
            v = v * coeff;
            return coeff;
    }

    void update_position_internal(NT&){}

    template <class bfunc, class NonLinearOracle>
    std::tuple<NT, Point, int> curve_intersect(
      NT t_prev,
      NT t0,
      NT eta,
      std::vector<Point> &coeffs,
      bfunc phi,
      bfunc grad_phi,
      NonLinearOracle &intersection_oracle,
      int ignore_facet=-1)
    {
        return intersection_oracle.apply(t_prev, t0, eta, A, b, *this,
                                         coeffs, phi, grad_phi, ignore_facet);
    }
};

#endif