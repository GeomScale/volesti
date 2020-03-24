// volesti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BALL_H
#define BALL_H

// ball type
template <typename Point>
struct Ball{
public:
    typedef Point BallPoint;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    Ball() {}

    Ball(Point cc, NT RR) : c(cc),	 R(RR) {}

    Point center() const {
        return c;
    }

    NT squared_radius() const {
        return R;
    }

    NT radius() const {
        return std::sqrt(R);
    }

    int dimension() const {
        return c.dimension();
    }

    int is_in(Point &p) {
        if (p.squared_length() <= R)
            return -1;
        else return 0;
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v) {

        NT vrc(0), v2(0), rc2(0);

        vrc = v.dot(r);
        v2 = v.dot(v);
        rc2 = r.dot(r);

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        return std::pair<NT,NT> ((NT(-1)*vrc + disc_sqrt)/v2, (NT(-1)*vrc - disc_sqrt)/v2);
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, const VT &Ar, const VT &Av){
        return line_intersect(r, v);
    }


    std::pair<NT,NT> line_intersect(Point &r, Point &v, const VT &Ar, const VT &Av, NT &lambda_prev) {
        return line_intersect(r, v);
    }

    std::pair<NT,int> line_positive_intersect(Point &r, Point &v){
        return std::pair<NT,NT>(line_intersect(r, v).first, 0);
    }

    std::pair<NT,int> line_positive_intersect(Point &r, Point &v, const VT &Ar,
                                             const VT &Av){
        return line_positive_intersect(r, v);
    }

    std::pair<NT,int> line_positive_intersect(Point &r, Point &v, const VT &Ar,
                                             const VT &Av, NT &lambda_prev){
        return line_positive_intersect(r, v);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r, const unsigned int &rand_coord) {

        NT vrc = r[rand_coord];
        NT rc2(R);
        rc2 -=  r.dot(r);


        NT disc_sqrt = std::sqrt(std::pow(vrc,2) + rc2);
        return std::pair<NT,NT> (NT(-1)*vrc + disc_sqrt, NT(-1)*vrc - disc_sqrt);

    }

    std::pair<NT,NT> line_intersect_coord(Point &r, const unsigned int &rand_coord,
                                          const VT &lamdas) {
        return line_intersect_coord(r, rand_coord);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const Point &r_prev,
                                          const unsigned int &rand_coord,
                                          const unsigned int &rand_coord_prev,
                                          const VT &lamdas) {
        return line_intersect_coord(r, rand_coord);
    }

    int num_of_hyperplanes() {
        return 0;
    }

    void compute_reflection (Point &v, const Point &p, const int &facet) {

        Point s = p;
        s *= (1.0 / std::sqrt(s.squared_length()));
        s *= (-2.0 * v.dot(s));
        v += s;
    }

private:
    Point  c; //center
    NT     R; //SQUARED radius !!!
};


#endif
