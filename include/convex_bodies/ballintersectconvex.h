// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BALLINTERSECTCONVEX_H
#define BALLINTERSECTCONVEX_H


template <typename Polytope, typename CBall>
class BallIntersectPolytope {
private:
    Polytope    P;
    CBall B;
public:
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;
    typedef typename CBall::NT NT;
    typedef typename CBall::PointType PointType;

    BallIntersectPolytope() {}

    BallIntersectPolytope(Polytope const& PP, CBall &BB) : P(PP), B(BB) {};

    Polytope first() const { return P; }
    CBall second() const { return B; }

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
        if (B.is_in(p)==-1)
            return P.is_in(p);
        return 0;
    }

    int num_of_hyperplanes() const {
        return P.num_of_hyperplanes();
    }

    unsigned int dimension() const {
        return P.dimension();
    }

    NT radius() const {
        return B.radius();
    }

    std::pair<NT,NT> line_intersect(PointType const& r, PointType const& v) const
    {

        std::pair <NT, NT> polypair = P.line_intersect(r, v);
        std::pair <NT, NT> ballpair = B.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT &Ar,
                                    VT &Av) const
    {
        std::pair <NT, NT> polypair = P.line_intersect(r, v, Ar, Av);
        std::pair <NT, NT> ballpair = B.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT &Ar,
                                    VT &Av,
                                    NT &lambda_prev) const
    {
        std::pair <NT, NT> polypair = P.line_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, NT> ballpair = B.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,int> line_positive_intersect(PointType const& r,
                                              PointType const& v,
                                              VT &Ar,
                                              VT &Av) const
    {

        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av);
        std::pair <NT, int> ball_lambda = B.line_positive_intersect(r, v);
        int facet = P.num_of_hyperplanes();

        if (polypair.first < ball_lambda.first ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ball_lambda.first), facet);
    }


    std::pair<NT,int> line_positive_intersect(PointType const& r,
                                              PointType const& v,
                                              VT &Ar,
                                              VT &Av,
                                              NT &lambda_prev) const
    {

        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, int> ball_lambda = B.line_positive_intersect(r, v);
        int facet = P.num_of_hyperplanes();

        if (polypair.first < ball_lambda.first ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ball_lambda.first), facet);
    }

    //First coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(PointType const& r,
                                          unsigned int const& rand_coord,
                                          VT &lamdas) const
    {

        std::pair <NT, NT> polypair = P.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> ballpair = B.line_intersect_coord(r, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
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
        std::pair <NT, NT> ballpair = B.line_intersect_coord(r, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,NT> query_dual(PointType const& p,
                                unsigned int const& rand_coord)
    {
        std::pair <NT, NT> polypair = P.query_dual(p, rand_coord);
        std::pair <NT, NT> ballpair = B.line_intersect_coord(p, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    void compute_reflection (PointType& v, PointType const& p, int &facet) const
    {

        if (facet == P.num_of_hyperplanes()) {
            B.compute_reflection(v, p, facet);
        } else {
            P.compute_reflection(v, p, facet);
        }

    }

};


/* EXPERIMENTAL
template <class T1 , class T2>
class PolytopeIntersectEllipsoid {
private:
    T1 P;
    T2 E;
    typedef typename T2::K 	K;
public:
    PolytopeIntersectEllipsoid(T1 &Pin, T2 &Ein) : P(Pin), E(Ein) {};

    T1 first() { return P; }
    T2 second() { return E; }

    int is_in(Point p){
        //std::cout << "calling is in"<<std::endl;
        if(P.is_in(p)==-1)
            return E.is_in(p);
        return 0;
    }

    int num_of_hyperplanes(){
        return P.num_of_hyperplanes();
    }

    int dimension(){
        return P.dimension();
    }

    std::pair<Point,Point> line_intersect(Point r,
                                          Point v){

        std::pair<Point,Point> polypair = P.line_intersect(r,v);
        std::pair<Point,Point> returnpair;
        std::pair<Point,Point> ellpair;
        bool ellinter=false;

        //check the first intersection point if it is inside ball
        if(E.is_in(polypair.first)){
            returnpair.first = polypair.first;
        }else{
            ellinter=true;
            //compute the intersection with ball
            ellpair = E.line_intersect(r,v);
            returnpair.first = ellpair.first;
        }
        //check the second intersection point
        if(E.is_in(polypair.second)){
            returnpair.second = polypair.second;
        }else{
            if(ellinter) //if the intersection with ball is already computed
                returnpair.second = ellpair.second;
            else returnpair.second = (E.line_intersect(r,v)).second;
        }
        return returnpair;
    }

    std::pair<K,K> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init
                                          ){

        std::pair<K,K> polypair = P.line_intersect_coord(r,r_prev,rand_coord,rand_coord_prev,lamdas,init);
        std::pair<K,K> ellpair = E.line_intersect_coord(r,rand_coord);
        return std::pair<K,K> (std::min(polypair.first,ellpair.first),
                                 std::max(polypair.second,ellpair.second));
    }

};


template <class T1 , class T2 >
class BallPolyIntersectEll {
private:
    T1 BP;
    T2 E;
    typedef typename T2::K 	K;
public:
    BallPolyIntersectEll(T1 &BPin, T2 &Ein) : BP(BPin), E(Ein) {};

    T1 first() { return BP; }
    T2 second() { return E; }

    int is_in(Point p){
        //std::cout << "calling is in"<<std::endl;
        if(BP.is_in(p)==-1)
            return E.is_in(p);
        return 0;
    }

    int num_of_hyperplanes(){
        return BP.num_of_hyperplanes();
    }

    int dimension(){
        return BP.dimension();
    }

    std::pair<Point,Point> line_intersect(Point r,
                                          Point v){

        std::pair<Point,Point> Bpolypair = BP.line_intersect(r,v);
        std::pair<Point,Point> returnpair;
        std::pair<Point,Point> ellpair;
        bool ellinter=false;

        //check the first intersection point if it is inside ball
        if(E.is_in(Bpolypair.first)){
            //std::cout<<"inside ball 1, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.first<<std::endl;
            returnpair.first = Bpolypair.first;
        }else{
            //std::cout<<"outside ball 1, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.first<<std::endl;
            ellinter=true;
            //compute the intersection with ball
            ellpair = E.line_intersect(r,v);
            returnpair.first = ellpair.first;
            //std::cout<<returnpair.first<<std::endl;
        }
        //check the second intersection point
        if(E.is_in(Bpolypair.second)){
            //std::cout<<"inside ball 2, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.second<<std::endl;
            returnpair.second = Bpolypair.second;
        }else{
            //std::cout<<"outside ball 2, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.second<<std::endl;
            if(ellinter) //if the intersection with ball is already computed
                returnpair.second = ellpair.second;
            else returnpair.second = (E.line_intersect(r,v)).second;
            //std::cout<<returnpair.second<<std::endl;
        }
        return returnpair;
    }

    std::pair<K,K> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init
                                          ){

        std::pair<K,K> Bpolypair = BP.line_intersect_coord(r,r_prev,rand_coord,rand_coord_prev,lamdas,init);
        std::pair<K,K> ellpair = E.line_intersect_coord(r,rand_coord);
        return std::pair<K,K> (std::min(Bpolypair.first,ellpair.first),
                                 std::max(Bpolypair.second,ellpair.second));
    }



};*/

#endif
