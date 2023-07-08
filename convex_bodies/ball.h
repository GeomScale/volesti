// volesti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018-19 programs.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BALL_H
#define BALL_H

#include <Eigen/Eigen>

/// This class represents a ball parameterized by a point type
/// \tparam Point Point Type
template <typename Point>
class Ball{
public:
    typedef Point PointType;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    Ball() {}

    Ball(Point cc, NT RR) : c(cc),	 R(RR) {}

    std::pair<Point,NT> InnerBall() const
    {
        return std::pair<Point,NT>(c, R);
    }

    Point center() const
    {
        return c;
    }

    NT squared_radius() const
    {
        return R;
    }

    NT radius() const
    {
        return std::sqrt(R);
    }

    int dimension() const
    {
        return c.dimension();
    }

    int is_in(Point const& p) const
    {
        if (p.squared_length() <= R)
            return -1;
        else return 0;
    }

    std::pair<NT,NT> line_intersect(Point const& r, Point const& v) const
    {

        NT vrc(0), v2(0), rc2(0);

        vrc = v.dot(r);
        v2 = v.dot(v);
        rc2 = r.dot(r);

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        return std::pair<NT,NT> ((NT(-1)*vrc + disc_sqrt)/v2,
                                 (NT(-1)*vrc - disc_sqrt)/v2);
    }

    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    const VT &Ar,
                                    const VT &Av) const
    {
        return line_intersect(r, v);
    }


    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    const VT &Ar,
                                    const VT &Av,
                                    NT &lambda_prev) const
    {
        return line_intersect(r, v);
    }

    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v) const
    {
        return std::pair<NT,NT>(line_intersect(r, v).first, 0);
    }

    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              const VT &Ar,
                                              const VT &Av) const
    {
        return line_positive_intersect(r, v);
    }

    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              const VT &Ar,
                                              const VT &Av,
                                              NT &lambda_prev) const
    {
        return line_positive_intersect(r, v);
    }

    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord) const
    {

        NT vrc = r[rand_coord];
        NT rc2(R);
        rc2 -=  r.dot(r);


        NT disc_sqrt = std::sqrt(std::pow(vrc,2) + rc2);
        return std::pair<NT,NT> (NT(-1)*vrc + disc_sqrt, NT(-1)*vrc - disc_sqrt);

    }

    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          const VT &lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }

    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          const VT &lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }

    int num_of_hyperplanes() const
    {
        return 0;
    }

    void compute_reflection (Point& v, Point const& p) const
    {
        Point s = p;
        s *= (1.0 / std::sqrt(s.squared_length()));
        s *= (-2.0 * v.dot(s));
        v += s;
    }

    template <typename update_parameters>
    void compute_reflection (Point &v, Point const& p, update_parameters &params) const {

        params.ball_inner_norm = p.length();
        params.inner_vi_ak = v.dot(p) / params.ball_inner_norm;
        v += (p * (-2.0 * params.inner_vi_ak * (1.0 / params.ball_inner_norm)));
    }

private:
    Point  c; //center
    NT     R; //SQUARED radius !!!
};


#endif
