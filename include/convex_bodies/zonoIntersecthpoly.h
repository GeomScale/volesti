// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

#ifndef ZONOINTERSECTHPOLY_H
#define ZONOINTERSECTHPOLY_H


template <typename Zonotope, typename HPolytope>
class ZonoIntersectHPoly {
private:
    Zonotope    Z;
    HPolytope HP;
public:
    typedef typename HPolytope::NT NT;
    typedef typename HPolytope::VT VT;
    typedef typename HPolytope::MT MT;
    typedef typename Zonotope::PolytopePoint Point;

    ZonoIntersectHPoly() {}

    ZonoIntersectHPoly(Zonotope &Z1, HPolytope &HP1) : Z(Z1), HP(HP1) {};

    Zonotope first() const { return Z; }
    HPolytope second() const { return HP; }

    int is_in(const Point &p){
        if(HP.is_in(p)==-1)
            return Z.is_in(p);
        return 0;
    }

    int num_of_hyperplanes() const {
        return HP.num_of_hyperplanes();
    }

    unsigned int dimension() const {
        return HP.dimension();
    }

    unsigned int num_of_generators() const {
        return Z.num_of_generators();
    }

    NT radius() const {
        return HP.radius();
    }

    MT get_mat() const {
        return HP.get_mat();
    }

    MT get_vec() const {
        return HP.get_vec();
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v) {

        std::pair <NT, NT> polypair = HP.line_intersect(r, v);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, VT &Ar,
            VT &Av) {
        std::pair <NT, NT> polypair = HP.line_intersect(r, v, Ar, Av);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, VT &Ar,
            VT &Av, NT &lambda_prev) {
        std::pair <NT, NT> polypair = HP.line_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    //First coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(Point &r,const unsigned int &rand_coord,
                                          VT &lamdas) {

        std::pair <NT, NT> polypair = HP.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }


    std::pair<NT, int> line_positive_intersect(Point &r, Point &v, VT &Ar,
                                               VT &Av) {

        std::pair <NT, int> polypair = HP.line_positive_intersect(r, v, Ar, Av);
        std::pair <NT, int> zonopair  = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes()+1;

        if (polypair.first < zonopair.first ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first), facet);
    }

    std::pair<NT, int> line_positive_intersect(Point &r, Point &v,  VT &Ar,
                                               VT &Av, NT &lambda_prev) {
        std::pair <NT, int> polypair = HP.line_positive_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, int> zonopair  = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes()+1;

        if (polypair.first < zonopair.first ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first), facet);
    }

    //Not the first coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const Point &r_prev,
                                          const unsigned int &rand_coord,
                                          const unsigned int &rand_coord_prev,
                                          VT &lamdas) {

        std::pair <NT, NT> polypair = HP.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> query_dual(Point &p, const unsigned int &rand_coord) {
        std::pair <NT, NT> polypair = HP.query_dual(p, rand_coord);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(p, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    void compute_reflection (Point &v, const Point &p, const int &facet) {

        if (facet == (HP.num_of_hyperplanes()+1)) {
            Z.compute_reflection(v, p, facet);
        } else {
            HP.compute_reflection(v, p, facet);
        }

    }

};

#endif
