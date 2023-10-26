// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2023 Vissarion Fisikopoulos
// Copyright (c) 2018-2023 Apostolos Chalkis
// Copyright (c) 2029-2023 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ELLIPSOIDINTERSECTCONVEX_H
#define ELLIPSOIDINTERSECTCONVEX_H

/// This class represents a polytope intersected with an ellipsoid
/// \tparam Polytope Polytope Type
/// \tparam CEllipsoid Ellipsoid Type
template <typename Polytope, typename CEllipsoid>
class EllipsoidIntersectPolytope {
private:
    Polytope    P;
    CEllipsoid E;
public:
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;
    typedef typename CEllipsoid::NT NT;
    typedef typename CEllipsoid::PointType PointType;

    EllipsoidIntersectPolytope() {}

    EllipsoidIntersectPolytope(Polytope& PP, CEllipsoid &EE) : P(PP), E(EE) {};

    Polytope first() const { return P; }
    CEllipsoid second() const { return E; }

    std::pair<PointType,NT> InnerBall() const
    {
        return P.InnerBall();
    }

    MT get_mat() const {
        return P.get_mat();
    }

    MT get_T() const {
        return P.get_mat();
    }

    MT get_vec() const {
        return P.get_vec();
    }

    int is_in(PointType const& p) const
    {
        if (P.is_in(p)==-1)
            return E.is_in(p);
        return 0;
    }

    int num_of_hyperplanes() const {
        return P.num_of_hyperplanes();
    }

    unsigned int dimension() const {
        return P.dimension();
    }

    NT radius() const {
        return E.radius();
    }

    std::pair<NT,NT> line_intersect(PointType const& r, PointType const& v) const
    {

        std::pair <NT, NT> polypair = P.line_intersect(r, v);
        std::pair <NT, NT> ellipsoidpair = E.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ellipsoidpair.first),
                                 std::max(polypair.second, ellipsoidpair.second));
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT &Ar,
                                    VT &Av) const
    {
        std::pair <NT, NT> polypair = P.line_intersect(r, v, Ar, Av);
        std::pair <NT, NT> ellipsoidpair = E.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ellipsoidpair.first),
                                 std::max(polypair.second, ellipsoidpair.second));
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT &Ar,
                                    VT &Av,
                                    NT &lambda_prev) const
    {
        std::pair <NT, NT> polypair = P.line_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, NT> ellipsoidpair = E.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ellipsoidpair.first),
                                 std::max(polypair.second, ellipsoidpair.second));
    }

    std::pair<NT,int> line_positive_intersect(PointType const& r,
                                              PointType const& v,
                                              VT &Ar,
                                              VT &Av)
    {
        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av);
        std::pair <NT, int> ellipsoid_lambda = E.line_positive_intersect(r, v);
        int facet = (polypair.first < ellipsoid_lambda.first) ? polypair.second : P.num_of_hyperplanes();

        return std::pair<NT, int>(std::min(polypair.first, ellipsoid_lambda.first), facet);
    }


    std::pair<NT,int> line_positive_intersect(PointType const& r,
                                              PointType const& v,
                                              VT &Ar,
                                              VT &Av,
                                              NT &lambda_prev)
    {
        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, int> ellipsoid_lambda = E.line_positive_intersect(r, v);
        int facet = (polypair.first < ellipsoid_lambda.first) ? polypair.second : P.num_of_hyperplanes();

        return std::pair<NT, int>(std::min(polypair.first, ellipsoid_lambda.first), facet);
    }

    //---------------------accelerated billiard---------------------//
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType const& r,
                                                     PointType const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters& params)
    {
        std::pair <NT, int> polypair = P.line_first_positive_intersect(r, v, Ar, Av, params);
        std::pair <NT, int> ellipsoid_lambda = E.line_positive_intersect(r, v);

        params.hit_ball = (polypair.first < ellipsoid_lambda.first) ? false : true;
        int facet = params.hit_ball ? P.num_of_hyperplanes() : polypair.second;
        params.facet_prev = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ellipsoid_lambda.first), facet);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               MT const& AA,
                                               update_parameters& params)
    {
        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av, lambda_prev, params);
        std::pair <NT, int> ellipsoid_lambda = E.line_positive_intersect(r, v);

        params.hit_ball = (polypair.first < ellipsoid_lambda.first) ? false : true;
        int facet = params.hit_ball ? P.num_of_hyperplanes() : polypair.second;
        params.facet_prev = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ellipsoid_lambda.first), facet);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters& params)
    {
        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av, lambda_prev, params);
        std::pair <NT, int> ellipsoid_lambda = E.line_positive_intersect(r, v);

        params.hit_ball = (polypair.first < ellipsoid_lambda.first) ? false : true;
        int facet = params.hit_ball ? P.num_of_hyperplanes() : polypair.second;
        params.facet_prev = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ellipsoid_lambda.first), facet);
    }
//-------------------------------------------------------------------------//

    //First coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(PointType const& r,
                                          unsigned int const& rand_coord,
                                          VT &lamdas) const
    {

        std::pair <NT, NT> polypair = P.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> ellipsoidpair = E.line_intersect_coord(r, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ellipsoidpair.first),
                                 std::max(polypair.second, ellipsoidpair.second));
    }

    //Not the first coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(PointType const& r,
                                          PointType const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          VT &lamdas) const
    {
        std::pair <NT, NT> polypair = P.line_intersect_coord(r, r_prev, rand_coord,
                                                             rand_coord_prev, lamdas);
        std::pair <NT, NT> ellipsoidpair = E.line_intersect_coord(r, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ellipsoidpair.first),
                                 std::max(polypair.second, ellipsoidpair.second));
    }

    std::pair<NT,NT> query_dual(PointType const& p,
                                unsigned int const& rand_coord)
    {
        std::pair <NT, NT> polypair = P.query_dual(p, rand_coord);
        std::pair <NT, NT> ellipsoidpair = E.line_intersect_coord(p, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ellipsoidpair.first),
                                 std::max(polypair.second, ellipsoidpair.second));
    }

    void compute_reflection (PointType& v, PointType const& p, int &facet)
    {

        if (facet == P.num_of_hyperplanes()) {
            E.compute_reflection(v, p);
        } else {
            P.compute_reflection(v, p, facet);
        }

    }

    template <typename update_parameters>
    void compute_reflection (PointType &v, PointType const& p, update_parameters &params)
    {
        if (params.hit_ball) {
            E.compute_reflection(v, p, params);
        } else {
            P.compute_reflection(v, p, params);
        }
    }

    void update_position_internal(NT&){}

    void resetFlags() {}

    std::pair<PointType, NT> ComputeInnerBall()
    {
        return P.ComputeInnerBall();
    }

};

#endif // ELLIPSOIDINTERSECTCONVEX_H
