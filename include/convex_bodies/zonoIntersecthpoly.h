// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018-19 programs.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#ifndef ZONOINTERSECTHPOLY_H
#define ZONOINTERSECTHPOLY_H

/// This class represents the intersection of a Zonotope with an H-polytope
/// \tparam Zonotope Zonotope Type
/// \tparam HPolytope HPolytope Type
template <typename Zonotope, typename HPolytope>
class ZonoIntersectHPoly {
private:
    Zonotope    Z;
    HPolytope HP;
public:
    typedef typename HPolytope::NT NT;
    typedef typename HPolytope::VT VT;
    typedef typename HPolytope::MT MT;
    typedef typename Zonotope::PointType PointType;

    ZonoIntersectHPoly()
    {}

    ZonoIntersectHPoly(Zonotope &Z1, HPolytope &HP1)
        : Z(Z1)
        , HP(HP1)
    {}

    Zonotope first() const { return Z; }
    HPolytope second() const { return HP; }

    int is_in(const PointType &p) const {
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

    MT get_mat() const {
        return HP.get_mat();
    }

    MT get_T() const {
        return Z.get_mat();
    }

    MT get_vec() const {
        return HP.get_vec();
    }

    std::pair<PointType,NT> InnerBall() const
    {
        return Z.InnerBall();
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v) const
    {

        std::pair <NT, NT> polypair = HP.line_intersect(r, v);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT& Ar,
                                    VT& Av) const
    {
        std::pair <NT, NT> polypair = HP.line_intersect(r, v, Ar, Av);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT& Ar,
                                    VT& Av,
                                    NT const& lambda_prev) const
    {
        std::pair <NT, NT> polypair = HP.line_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, NT> zonopair = Z.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    //First coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(PointType const& r,
                                          unsigned int const& rand_coord,
                                          VT& lamdas) const
    {
        std::pair <NT, NT> polypair = HP.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(r, rand_coord, lamdas);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }


    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT& Ar,
                                               VT& Av) const
    {
        std::pair <NT, int> polypair = HP.line_positive_intersect(r, v, Ar, Av);
        std::pair <NT, int> zonopair  = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes()+1;

        if (polypair.first < zonopair.first ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first), facet);
    }

    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev) const
    {
        std::pair <NT, int> polypair = HP.line_positive_intersect(r, v, Ar, Av,
                                                                  lambda_prev);
        std::pair <NT, int> zonopair  = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes()+1;

        if (polypair.first < zonopair.first ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first),
                                  facet);
    }

    //---------------------accelarated billiard---------------------//
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType const& r,
                                                     PointType const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters& params) const
    {
        std::pair <NT, int> polypair = HP.line_first_positive_intersect(r, v, Ar, Av, params);
        std::pair <NT, int> zonopair = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes();
        params.facet_prev = polypair.second;

        if (polypair.first < zonopair.first ) {
            facet = polypair.second;
            params.hit_ball = false;
        } else {
            params.hit_ball = true;
        }

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first), facet);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               MT const& AA,
                                               update_parameters& params) const
    {
        std::pair <NT, int> polypair = HP.line_positive_intersect(r, v, Ar, Av, lambda_prev, params);
        std::pair <NT, int> zonopair = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes();
        params.facet_prev = polypair.second;

        if (polypair.first < zonopair.first ) {
            facet = polypair.second;
            params.hit_ball = false;
        } else {
            params.hit_ball = true;
        }

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first), facet);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters& params) const
    {
        std::pair <NT, int> polypair = HP.line_positive_intersect(r, v, Ar, Av, lambda_prev, params);
        std::pair <NT, int> zonopair = Z.line_positive_intersect(r, v, Ar, Av);
        int facet = HP.num_of_hyperplanes();
        params.facet_prev = polypair.second;

        if (polypair.first < zonopair.first ) {
            facet = polypair.second;
            params.hit_ball = false;
        } else {
            params.hit_ball = true;
        }

        return std::pair<NT, int>(std::min(polypair.first, zonopair.first), facet);
    }
//-------------------------------------------------------------------------//

    //Not the first coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(PointType const& r,
                                          PointType const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          VT& lamdas) const
    {

        std::pair <NT, NT> polypair = HP.line_intersect_coord(r, r_prev,
                                                              rand_coord,
                                                              rand_coord_prev,
                                                              lamdas);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(r, rand_coord,
                                                             lamdas);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    std::pair<NT,NT> query_dual(PointType const& p,
                                unsigned int const& rand_coord) const
    {
        std::pair <NT, NT> polypair = HP.query_dual(p, rand_coord);
        std::pair <NT, NT> zonopair = Z.line_intersect_coord(p, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, zonopair.first),
                                 std::max(polypair.second, zonopair.second));
    }

    void compute_reflection (PointType &v,
                             PointType const& p,
                             int const& facet) const
    {
        if (facet == (HP.num_of_hyperplanes()+1)) {
            Z.compute_reflection(v, p, facet);
        } else {
            HP.compute_reflection(v, p, facet);
        }

    }

    template <typename update_parameters>
    void compute_reflection (PointType &v, const PointType &p, update_parameters const& params) const {

        if (params.hit_ball) {
            Z.compute_reflection(v, p, params);
        } else {
            HP.compute_reflection(v, p, params.facet_prev);
        }

    }

};

#endif
