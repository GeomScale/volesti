// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HPOLYTOPE_H
#define HPOLYTOPE_H

#include <limits>
#include <iostream>
#include <Eigen/Eigen>
#include "lp_oracles/solve_lp.h"
#include "nlp_oracles/nlp_hpolyoracles.h"

#define MAX_NR_TRIES 10000

//min and max values for the Hit and Run functions


// H-polytope class
template <typename Point>
class HPolytope{
public:
    typedef Point PointType;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    //using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    //typedef RowMatrixXd MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    MT A; //matrix A
    VT b; // vector b, s.t.: Ax<=b
    unsigned int            _d; //dimension
    std::pair<Point,NT> _inner_ball;
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();

public:

    HPolytope() {}

    std::pair<Point,NT> InnerBall() const
    {
        return _inner_ball;
    }

    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point,NT> ComputeInnerBall()
    {
        _inner_ball = ComputeChebychevBall<NT, Point>(A, b);
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


    // return the vector b
    VT get_vec() const
    {
        return b;
    }


    // change the matrix A
    void set_mat(MT const& A2)
    {
        A = A2;
    }


    // change the vector b
    void set_vec(VT const& b2)
    {
        b = b2;
    }

    Point get_mean_of_vertices() const
    {
        return Point(_d);
    }

    NT get_max_vert_norm() const
    {
        return 0.0;
    }

    void init(unsigned int const& dim, MT const& _A, VT const& _b)
    {
        _d = dim;
        A = _A;
        b = _b;
    }

    //define matrix A and vector b, s.t. Ax<=b and the dimension
    void init(std::vector<std::vector<NT> > const& Pin)
    {
        _d = Pin[0][1] - 1;
        A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                A(i - 1, j - 1) = -Pin[i][j];
            }
        }
    }


    // print polytope in input format
    void print() {
        std::cout << " " << A.rows() << " " << _d << " float" << std::endl;
        for (unsigned int i = 0; i < A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                std::cout << A(i, j) << " ";
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
    int is_in(Point const& p) const
    {
        int m = A.rows();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            //Check if corresponding hyperplane is violated
            if (*b_data - A.row(i) * p.getCoefficients() < NT(0))
                return 0;

            b_data++;
        }
        return -1;
    }


    template <class bfunc>
    std::tuple<NT, Point, int> curve_intersect_ipopt(NT t_prev, NT t0, std::vector<Point> &coeffs, bfunc phi, bfunc grad_phi) {
      return curve_intersect_hpoly_ipopt_helper<MT, VT, Point, NT, bfunc>(t_prev, t0, A, b, coeffs, phi, grad_phi);
    }

    // Compute intersection of H-polytope P := Ax <= b
    // with curve p(t) = sum a_j phi_j(t) where phi_j are basis
    // functions (e.g. polynomials)
    // Uses Newton-Raphson to solve the transcendental equation
    template <class bfunc>
    std::tuple<NT, Point, int> curve_intersect_newton_raphson(NT t_prev, NT t0, std::vector<Point> &coeffs, bfunc phi, bfunc grad_phi) {

      // Keep results in a vector (in case of multiple roots)
      // The problem has O(m * len(coeffs)) solutions if phi's are polys
      // due to the Fundamental Theorem of Algebra
      // Some roots may be common for more than one hyperplanes
      // The equations may have complex roots as well but they do not
      // interest us (we don't find them)
      // std::vector<std::tuple<NT, Point, int>> results;

      std::tuple<NT, Point, int> result = std::make_tuple(NT(maxNT), Point(_d), -1);

      // Root
      NT t = t_prev;

      // Helper variables for Newton-Raphson
      NT dot_u, num, den, den_tmp;

      // Regularization for NR (e.g. at critical points where grad = 0)
      NT reg = (NT) 1e-7;
      VT u, Z;
      int m = num_of_hyperplanes();

      // Keeps constants A_i^T C_j
      Z.resize(coeffs.size());

      // Helper vector (lies on m-th hyperplane)
      u.resize(_d);

      // Iterate over all hyperplanes
      for (int i = 0; i < m; i++) {

        // Calculate constants
        start_iter: t_prev = t0 + reg;


        for (unsigned int j = 0; j < coeffs.size(); j++) {
          Z(j) = A.row(i) * coeffs[j].getCoefficients();
        }

        // Find point on m-th hyperplane
        for (unsigned int j = 0; j < _d; j++) u(j, 0) = 0;

        // If b[i] = 0 then point is (0, 0, ..., 0)
        if (!(b(i) == 0)) {
          for (unsigned int j = 0; j < _d; j++) {
            // Else A(i) must have a non-zero entry
            // Find it and set the coefficient equal to A(i, j) / b(i)
            // Set the others to 0
            if (!(A(i, j) == 0)) {
              u(j, 0) = b(i) / A(i, j);
              break;
            }
          }
          // The point (0, 0, ..., A(i, j) / b(i), 0, 0, ... ) is on the m-th hyperplane
        }

        NT dot_u = (NT) (A.row(i) * u.col(0));

        for (int tries = 0; tries < MAX_NR_TRIES; tries++) {

          num = - dot_u;
          den = (NT) 0;

          // Calculate numerator f(t) and denominator f'(t)
          for (int j = 0; j < coeffs.size(); j++) {
            num += Z(j) * phi(t_prev, t0, j, coeffs.size());

            // Avoid ill-posed derivative (e.g. 0^{-1})
            if (j > 0) den += Z(j) * grad_phi(t_prev, t0, j, coeffs.size());
          }

          // Regularize denominator if near 0
          if (std::abs(den) < 10 * reg) den += reg;

          // Newton-Raphson Iteration t = t_old - f(t) / f'(t)
          t = t_prev - num / den;

          if (t < 0 && t_prev < 0) continue;

          if (std::abs(t - t_prev) < 1e-6 && t > 0) {
            // Add root (as t) and facet


            Point p = Point(coeffs[0].dimension());

            for (unsigned int j = 0; j < coeffs.size(); j++) {
              p += coeffs[j] * phi(t, t0, j, coeffs.size());
            }

            // TODO Keep largest positive root
            if (is_in(p) && t > 0 && t < std::get<0>(result)) result =  std::make_tuple(t, p, i);


          }

          t_prev = t;

        }

      }

      return result;

    }


    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point const& r, Point const& v) const
    {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
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
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT& Ar,
                                    VT& Av,
                                    bool pos = false) const
    {
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
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
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
    }

    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT& Ar,
                                    VT& Av,
                                    NT const& lambda_prev,
                                    bool pos = false) const
    {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
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
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
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


    //First coordinate ray intersecting convex polytope

    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          VT& lamdas) const
    {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
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
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          VT& lamdas) const
    {
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);

        int m = num_of_hyperplanes();

        lamdas.noalias() += A.col(rand_coord_prev)
                         * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
        NT* data = lamdas.data();

        for (int i = 0; i < m; i++) {
            NT a = A(i, rand_coord);

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
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    void linear_transformIt(MT const& T)
    {
        A = A * T;
    }


    // shift polytope by a point c

    void shift(const VT &c)
    {
        b -= A*c;
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
        NT row_norm;
        for (int i = 0; i < num_of_hyperplanes(); ++i)
        {
            row_norm = A.row(i).norm();
            A.row(i) = A.row(i) / row_norm;
            b(i) = b(i) / row_norm;
        }
    }

    void compute_reflection(Point& v, Point const&, int const& facet) const
    {
        v += -2 * v.dot(A.row(facet)) * A.row(facet);
    }

    void free_them_all() {}

};

#endif
