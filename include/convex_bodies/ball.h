// volesti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

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

    int is_in(Point &p) {
        if (p.squared_length() <= R)
            return -1;
        else return 0;
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v) {

        viterator vit=v.iter_begin(), cit=c.iter_begin(), rcit=r.iter_begin();
        NT vrc(0), v2(0), rc2(0);
        for( ; cit < c.iter_end() ; ++cit, ++rcit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        return std::pair<NT,NT> ((NT(-1)*vrc + disc_sqrt)/v2, (NT(-1)*vrc - disc_sqrt)/v2);
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, const std::vector<NT> &Ar,
            const std::vector<NT> &Av){
        return line_intersect(r, v);
    }


    std::pair<NT,NT> line_intersect(Point &r, Point &v, const std::vector<NT> &Ar,
            const std::vector<NT> &Av, NT &lambda_prev) {
        return line_intersect(r, v);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r, const unsigned int &rand_coord) {

        viterator rcit=r.iter_begin();
        NT vrc = *(rcit + rand_coord);
        NT rc2(R);
        for( ; rcit < r.iter_end() ; ++rcit){
            rc2 -= *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) + rc2);
        return std::pair<NT,NT> (NT(-1)*vrc + disc_sqrt, NT(-1)*vrc - disc_sqrt);

    }

    std::pair<NT,NT> line_intersect_coord(Point &r, const unsigned int &rand_coord,
                                          const std::vector<NT> &lamdas) {
        return line_intersect_coord(r, rand_coord);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const Point &r_prev,
                                          const unsigned int &rand_coord,
                                          const unsigned int &rand_coord_prev,
                                          const std::vector<NT> &lamdas) {
        return line_intersect_coord(r, rand_coord);
    }

    int num_of_hyperplanes() {
        return 0;
    }

private:
    Point  c; //center
    NT     R; //SQUARED radius !!!
};


#endif
