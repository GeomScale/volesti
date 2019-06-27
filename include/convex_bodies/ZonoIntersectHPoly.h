// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef ZONOINTERSECTHPOLY_H
#define ZONOINTERSECTHPOLY_H


template <class Zonotope, class HPolytope>
class ZonoIntersectHPoly {
private:
    Zonotope    Z;
    HPolytope HP;
public:
    typedef typename HPolytope::NT NT;
    typedef typename HPolytope::PolytopePoint Point;

    ZonoIntersectHPoly() {}

    ZonoIntersectHPoly(Zonotope &Z1, HPolytope &HP1) : Z(Z1), HP(HP1) {};

    Zonotope first() { return Z; }
    HPolytope second() { return HP; }

    int is_in(Point p){
        if(HP.is_in(p)==-1)
            return Z.is_in(p);
        return 0;
    }

    int num_of_hyperplanes(){
        return HP.num_of_hyperplanes();
    }

    unsigned int dimension(){
        return HP.dimension();
    }

    unsigned int num_of_generators(){
        return Z.num_of_generators();
    }

    std::pair<NT,NT> line_intersect(Point r, Point v) {

        std::pair <NT, NT> polypair = HP.line_intersect(r, v);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        std::pair <NT, NT> polypair = HP.line_intersect(r, v, Av, Ar);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        std::pair <NT, NT> polypair = HP.line_intersect(r, v, Av, Ar, lambda_prev);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    //First coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        std::pair <NT, NT> polypair = HP.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        return std::pair<NT, int> (0.0, 0);
    }

    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        return std::pair<NT, int> (0.0, 0);
    }

    //Not the first coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        std::pair <NT, NT> polypair = HP.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> query_dual(Point &p, unsigned int rand_coord) {
        std::pair <NT, NT> polypair = HP.query_dual(p, rand_coord);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(p, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    void compute_reflection (Point &v, Point &p, int &facet) {}

};

#endif
