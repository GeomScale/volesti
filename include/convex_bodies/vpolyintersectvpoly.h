// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VPOLYINTERSECTVPOLY_H
#define VPOLYINTERSECTVPOLY_H

template <class VPolytope>
class IntersectionOfVpoly {
public:
    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::PolytopePoint PolytopePoint;
    typedef PolytopePoint Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    VPolytope P1;
    VPolytope P2;


    IntersectionOfVpoly(VPolytope &P, VPolytope &Q) : P1(P), P2(Q) {};

    VPolytope first() { return P1; }
    VPolytope second() { return P2; }

    int is_in(Point p){
        if(P1.is_in(p)==-1)
            return P2.is_in(p);
        return 0;
    }

    int num_of_hyperplanes(){
        return 0;
    }

    unsigned int dimension() {
        return P1.dimension();
    }

    int num_of_vertices() {
        return P1.num_of_vertices() + P2.num_of_vertices();
    }

    unsigned int upper_bound_of_hyperplanes() {
        return P1.upper_bound_of_hyperplanes() + P2.upper_bound_of_hyperplanes() ;
    }

    std::pair<Point,NT> ComputeInnerBall() {

        unsigned int d = dimension();
        MT V(0, d);
        MT V1 = P1.get_mat();
        MT V2 = P2.get_mat();
        VT itervec(d);
        std::vector<NT> temp_p(d, 0.0);
        typename std::vector<NT>::iterator pit;
        int j;
        for (int i = 0; i < V1.rows(); ++i) {
            pit = temp_p.begin();
            j = 0;
            for (; pit!=temp_p.end(); ++pit, ++j) {
                *pit = V1(i,j);
                itervec(j) = V1(i,j);
            }
            Point p(d, temp_p.begin(), temp_p.end());
            if (P2.is_in(p) == -1) {
                V.conservativeResize(V.rows() + 1, V.cols());
                V.row(V.rows()-1) = itervec;
            }
        }

        for (int i = 0; i < V2.rows(); ++i) {
            pit = temp_p.begin();
            j = 0;
            for (; pit!=temp_p.end(); ++pit, ++j) {
                *pit = V2(i,j);
                itervec(j) = V2(i,j);
            }
            Point p(d, temp_p.begin(), temp_p.end());
            if (P1.is_in(p) == -1) {
                V.conservativeResize(V.rows() + 1, V.cols());
                V.row(V.rows()-1) = itervec;
            }
        }

        VPolytope Q;
        Q.init(d, V, itervec);
        return Q.ComputeInnerBall();
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r,
                                    Point v) {

        std::pair <NT, NT> P1pair = P1.line_intersect(r, v);
        std::pair <NT, NT> P2pair = P2.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));

    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {
        std::pair <NT, NT> P1pair = P1.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> P2pair = P2.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {
        std::pair <NT, NT> P1pair = P1.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        std::pair <NT, NT> P2pair = P2.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        return std::pair<NT, NT>(std::min(P1pair.first, P2pair.first),
                                 std::max(P1pair.second, P2pair.second));
    }


    // shift polytope by a point c
    void shift(VT c) {
        P1.shift(c);
        P2.shift(c);
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(MT T) {
        P1.linear_transformIt(T);
        P2.linear_transformIt(T);
    }

    std::vector<NT> get_dists(NT radius) {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    template <class PointList>
    bool get_points_for_rounding (PointList &randPoints) {
        if (num_of_vertices()>40*dimension()) {
            return false;
        }
        if(!P1.get_points_for_rounding(randPoints)) {
            return false;
        }
        if(!P2.get_points_for_rounding(randPoints)) {
            return false;
        }

        return true;
    }

};


#endif
