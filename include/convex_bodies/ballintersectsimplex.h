// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018-19 programs.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BALLINTERSECTSIMPLEX_H
#define BALLINTERSECTSIMPLEX_H

#include <limits>
#include <iostream>
#include <cmath>
#include <Eigen/Eigen>



//min and max values for the Hit and Run functions
// H-polytope class
template <typename NT, typename VT, typename MT>
class UnitBallIntersectSimplex {
public:
    //typedef Point                                             PointType;
    //typedef typename Point::FT                                NT;
    //typedef typename std::vector<NT>::iterator                viterator;
    //using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    //typedef RowMatrixXd MT;
    //typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    //typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;

private:
    unsigned int         _d; //dimension
    MT                   A; //matrix A
    VT                   b; // vector b, s.t.: Ax<=b
    MT                   V; // vertices of the simplex
    VT                   x0; // center of the unit sphere
    VT                   Vnorms;

public:
    //TODO: the default implementation of the Big3 should be ok. Recheck.
    UnitBallIntersectSimplex() {}

    UnitBallIntersectSimplex(unsigned d_, MT const& A_, VT const& b_, MT _V, VT _x0, VT _Vnorms) :
        _d{d_}, A{A_}, b{b_}, V{_V}, x0{_x0}, Vnorms{_Vnorms}
    {
    }

    // Copy constructor
    UnitBallIntersectSimplex(HPolytope<Point> const& p) :
            _d{p._d}, A{p.A}, b{p.b}, V{p.V}, x0{p.x0}, Vnorms{p.Vnorms}
    {
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

    void set_vertices(MT const& _V)
    {
        V = _V;
    }

    void set_center(VT const& _x0)
    {
        x0 = _x0;
    }

    // print polytope in input format
    void print() {
        std::cout << " " << A.rows() << " " << _d << " double" << std::endl;
        for (unsigned int i = 0; i < A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                std::cout << A(i, j) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
        }
    }


    int is_in(VT const& p, NT tol=NT(0)) const
    {
        int m = A.rows();
        VT temp = A * p - b;
        //VT temp2(_d);
        const NT* Ax_b_data = temp.data();
        for (int i = 0; i < m; i++) {
            //Check if corresponding hyperplane is violated
            if ((*Ax_b_data) < NT(-tol))
                return 0;

            Ax_b_data++;
        }

        VT v = p - x0;
        std::pair<NT,NT> pair_root = line_intersect(p, v);
        VT r = p + (pair_root.second * v);

        int n = V.cols();

        for (int i = 0; i < n; i++)
        {
            v = V.col(i) - r;
            r -= x0;

            NT a = v.dot(v);
            NT b = NT(2) * (r.dot(v));
            NT g = r.dot(r) - NT(1);

            NT D = b*b - NT(4) * a * g;

            if (D < NT(0))
            {
                return -1;
            }

            NT tmin = (-b - sqrt(D)) / (NT(2)*a);
            NT tmax = (-b + sqrt(D)) / (NT(2)*a);

            if (tmin < 1 && tmin > 0){
                continue;
            }
            else if (tmax < 1 && tmax > 0)
            {
                continue;
            }
            else
            {
                return -1;
            }
        }
        return 0;
    }

    // here the center x0 is always the origin
    int is_in_optimized(VT const& p, VT& Ar, VT& Av, NT &lambda_prev, NT tol=NT(0)) const
    {
        int m = A.rows();
        NT min_plus = std::numeric_limits<NT>::max();
        Ar.noalias() = cos(lambda_prev)*Ar + sin(lambda_prev)*Av;
        //VT b_Ax = b - Ar;
        //VT temp = A * p - b;
        //VT temp2(_d);
        NT* Ar_data = Ar.data();
        //NT* b_Ax_data = b_Ax.data();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            //Check if corresponding hyperplane is violated
            if ((*Ar_data) - (*b_data)< NT(-tol))
                return 0;

            NT lamda = ((*b_data) - (*Ar_data)) / (*Ar_data);
            if (lamda > NT(0) && lamda < min_plus)
            {
                min_plus = lamda;
            }
            Ar_data++;
            b_data++;
            b_Ax_data++;
        }

        VT r = p + (min_plus * v);

        int n = V.cols();

        for (int i = 0; i < n; i++)
        {
            NT r_v = r.dot(v);
            NT r_r = r.dot(r);

            NT a = Vnorms(i) - NT(2) * r_v + r_r;
            NT b = NT(2) * (r_v - r_r);
            NT g = r.dot(r) - NT(1);

            NT D = b*b - NT(4) * a * g;

            if (D < NT(0))
            {
                return -1;
            }

            NT tmin = (-b - sqrt(D)) / (NT(2)*a);
            NT tmax = (-b + sqrt(D)) / (NT(2)*a);

            if (tmin < 1 && tmin > 0){
                continue;
            }
            else if (tmax < 1 && tmax > 0)
            {
                continue;
            }
            else
            {
                return -1;
            }
        }
        return 0;
    }



    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(VT const& r, VT const& v) const
    {

        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom, sum_denom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes();


        sum_nom.noalias() = b - A * r;
        sum_denom.noalias() = A * v;

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
    std::pair<NT,NT> line_intersect(VT const& r,
                                    VT const& v,
                                    VT& Ar,
                                    VT& Av,
                                    bool pos = false) const
    {
        NT lamda = 0;
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        VT sum_nom;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() = A * r;
        sum_nom = b - Ar;
        Av.noalias() = A * v;


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

    std::pair<NT,NT> line_intersect(VT const& r,
                                    VT const& v,
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
        Av.noalias() = A * v;

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


   
    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> gc_intersect(Point const& r,
                                  Point const& v,
                                  VT& Ar,
                                  VT& Av) const
    {
        Ar.noalias() = A * r.getCoefficients();
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections(Ar, Av);
    }

    
    std::pair<NT,NT> gc_intersect(Point const& r,
                                    Point const& v,
                                    VT& Ar,
                                    VT& Av,
                                    NT const& lambda_prev) const
    {
        Ar.noalias() = cos(lambda_prev)*Ar + sin(lambda_prev)*Av;
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections(Ar, Av);
    }

    std::pair<NT,NT> gc_intersect_optimized(Point const& r,
                                            Point const& v,
                                            VT& Ar,
                                            VT& Av,
                                            NT const& lambda_prev) const
    {
        //Ar.noalias() = cos(lambda_prev)*Ar + sin(lambda_prev)*Av;
        Av.noalias() = A * v.getCoefficients();

        return compute_intersections(Ar, Av);
    }

    std::pair<NT,NT> compute_intersections(VT& Ar, VT& Av) const
    {
        NT D, C1, C2, eval, max_root = std::numeric_limits<NT>::lowest(), min_root = std::numeric_limits<NT>::max();
        NT min_plus  = std::numeric_limits<NT>::max();
        NT max_minus = std::numeric_limits<NT>::lowest();
        int m = num_of_hyperplanes();
        bool set_negative_root = false, set_positive_root = false, pos_D = false;

        NT* Av_data = Av.data();
        NT* Ar_data = Ar.data();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            D = (*Ar_data)*(*Ar_data) + (*Av_data)*(*Av_data) - (*b_data)*(*b_data);
            if (D > NT(0)) 
            {
                pos_D = true;

                C1 = asin((((*Av_data)*(*b_data)) + ((*Ar_data)*sqrt(D))) / ((*Ar_data)*(*Ar_data) + (*Av_data)*(*Av_data)));
                C2 = asin((((*Av_data)*(*b_data)) - ((*Ar_data)*sqrt(D))) / ((*Ar_data)*(*Ar_data) + (*Av_data)*(*Av_data)));

                eval = (*Ar_data) * cos(C1) + (*Av_data) * sin(C1) - (*b_data);
                if (!(eval > -NT(1e-05) && eval < NT(1e-05)))
                {
                    C1 = M_PI - C1;
                }
                if (C1 < min_plus && C1 > 0) {
                    min_plus = C1;
                    set_positive_root = true;
                }else if (C1 > max_minus && C1 < 0){
                    max_minus = C1;
                    set_negative_root = true;
                }
                if (C1 > max_root && C1 < NT(2)*M_PI)
                {
                    max_root = C1;
                }
                if (C1 < min_root && C1 > (-NT(2)*M_PI))
                {
                    min_root = C1;
                }

                eval = (*Ar_data) * cos(C2) + (*Av_data) * sin(C2) - (*b_data);
                if (!(eval > -NT(1e-05) && eval < NT(1e-05)))
                {
                    C2 = M_PI - C2;
                }
                if (C2 < min_plus && C2 > 0) {
                    min_plus = C2;
                    set_positive_root = true;
                }else if (C2 > max_minus && C2 < 0){
                    max_minus = C2;
                    set_negative_root = true;
                }
                if (C2 > max_root && C2 < NT(2)*M_PI)
                {
                    max_root = C2;
                }
                if (C2 < min_root && C2 > (-NT(2)*M_PI))
                {
                    min_root = C2;
                }
            }

            if (!set_negative_root)
            {
                if (pos_D)
                {
                    max_minus = max_root - NT(2) * M_PI;
                }
                else
                {
                    max_minus = NT(0);
                }
            }
            if (!set_positive_root)
            {
                if (pos_D)
                {
                    min_plus = min_root + NT(2) * M_PI;
                }
                else
                {
                    min_plus = NT(2) * M_PI;
                }
            }

            Av_data++;
            Ar_data++;
            b_data++;
        }
        return std::make_pair(min_plus, max_minus);
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

    void scale(const NT c)
    {
        b *= c;
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

};

#endif
