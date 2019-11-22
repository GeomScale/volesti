// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BALLINTERSECTCONVEX_H
#define BALLINTERSECTCONVEX_H

#include "solve_lp.h"

// ball type
template <class Point>
struct Ball{
public:
    typedef Point BallPoint;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;

    Ball() {}

    Ball(Point cc, NT RR) : c(cc),	 R(RR) {}

    Point center(){
        return c;
    }

    int dimension() {
        return c.dimension();
    }

    NT squared_radius(){
        return R;
    }

    NT radius(){
        return std::sqrt(R);
    }

    void set_radius(NT rad) {
        R = rad * rad;
    }

    int is_in(Point p){
        if (p.squared_length() <= R)
            return -1;
        else return 0;
    }

    int num_of_hyperplanes() {
        return 0;
    }

    std::pair<NT,NT> line_intersect(Point r, Point v){

        viterator rit=r.iter_begin();
        viterator vit=v.iter_begin();
        viterator cit=c.iter_begin();
        //Point rc = r;// - _c;
        viterator rcit=r.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        return std::pair<NT,NT> (lamda1,lamda2);
    }


    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av){

        viterator rit=r.iter_begin();
        viterator vit=v.iter_begin();
        viterator cit=c.iter_begin();
        //Point rc = r;// - _c;
        viterator rcit=r.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        return std::pair<NT,NT> (lamda1,lamda2);
    }


    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev){

        viterator rit=r.iter_begin();
        viterator vit=v.iter_begin();
        viterator cit=c.iter_begin();
        //Point rc = r;// - _c;
        viterator rcit=r.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        return std::pair<NT,NT> (lamda1,lamda2);
    }


    NT line_positive_intersect(Point r, Point v){

        viterator rit=r.iter_begin();
        viterator vit=v.iter_begin();
        viterator cit=c.iter_begin();

        viterator rcit=r.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        //NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        return lamda1;
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av){

        viterator rit=r.iter_begin();
        viterator vit=v.iter_begin();
        viterator cit=c.iter_begin();

        viterator rcit=r.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        //NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        return std::pair<NT,NT>(lamda1,0);
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev){

        viterator rit=r.iter_begin();
        viterator vit=v.iter_begin();
        viterator cit=c.iter_begin();

        viterator rcit=r.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        //NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        return std::pair<NT,NT>(lamda1,0);
    }


    std::pair<NT,NT> line_intersect_coord(Point r,
                                          int rand_coord){

        //Point rc = r;// - _c;
        viterator rcit=r.iter_begin();
        NT vrc = *(rcit + rand_coord);

        //NT v2 = NT(1);
        NT rc2(R);
        for( ; rcit < r.iter_end() ; ++rcit){
            rc2 -= *rcit * (*rcit);
        }

        //NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - _R));
        NT disc_sqrt = std::sqrt(std::pow(vrc,2) + rc2);// + _R);
        NT lamda1((NT(-1)*vrc + disc_sqrt));
        NT lamda2((NT(-1)*vrc - disc_sqrt));

        return std::pair<NT,NT> (lamda1,lamda2);

    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {
        return std::pair<NT,NT> (0.0, 0.0);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas){
        return std::pair<NT,NT> (0.0, 0.0);
    }

    void compute_reflection(Point &v, Point &p, int &facet) {

        Point s = (-1.0)*p;
        s = s * (1.0 / std::sqrt(s.squared_length()));
        s = ((-2.0 * v.dot(s)) * s);
        v = s + v;

    }

private:
    Point  c; //center
    NT     R; //SQUARED radius !!!
};


template <class Polytope, class CBall>
class BallIntersectPolytope {
private:
    Polytope    P;
    CBall B;
public:
    typedef typename CBall::NT NT;
    typedef typename CBall::BallPoint Point;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    typedef Polytope Hpoly;

    BallIntersectPolytope() {}

    BallIntersectPolytope(Polytope &PP, CBall &BB) : P(PP), B(BB) {};
    
    Polytope first() { return P; }
    CBall second() { return B; }

    int is_in(Point p){
        if(B.is_in(p)==-1)
            return P.is_in(p);
        return 0;
    }

    int num_of_hyperplanes(){
        return P.num_of_hyperplanes();
    }

    unsigned int dimension(){
        return P.dimension();
    }

    void set_vec(VT bb) {
        P.set_vec(bb);
    }

    VT get_vec() {
        return P.get_vec();
    }

    NT radius() {
        return B.radius();
    }

    void set_radius(NT rad) {
        B.set_radius(rad);
    }

    void comp_diam(NT &diam) {
        diam = 2.0*B.radius();
    }

    void add_facet(VT a, NT z0){
        P.add_facet(a, z0);
    }

    std::pair<Point,NT> ComputeInnerBall() {

        VT b = P.get_vec();
        int m = P.get_mat().rows(), d = P.get_mat().cols();
        NT rad = std::numeric_limits<NT>::max();
        std::cout<<"m = "<<m<<std::endl;
        for (int i = 0; i < m; ++i) {

            if (b(i) < rad) rad = b(i);
        }
        std::cout<<"CHeb radius of BP = "<<rad<<std::endl;
        return std::pair<Point, NT> (Point(d), rad);
        //return  ComputeChebychevBall2<Point>(P.get_mat(), P.get_vec(), radius());//ComputeChebychevBall<NT, Point>(A, b);
    }

    std::pair<NT,NT> line_intersect(Point r, Point v) {

        std::pair <NT, NT> polypair = P.line_intersect(r, v);
        std::pair <NT, NT> ballpair = B.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        std::pair <NT, NT> polypair = P.line_intersect(r, v, Ar, Av);
        std::pair <NT, NT> ballpair = B.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        std::pair <NT, NT> polypair = P.line_intersect(r, v, Ar, Av, lambda_prev);
        std::pair <NT, NT> ballpair = B.line_intersect(r, v);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av);
        NT ball_lambda = B.line_positive_intersect(r, v);
        int facet = P.num_of_hyperplanes();

        if (polypair.first < ball_lambda ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ball_lambda), facet);
    }


    std::pair<NT,int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av,
                                              NT &lambda_prev) {

        std::pair <NT, int> polypair = P.line_positive_intersect(r, v, Ar, Av, lambda_prev);
        NT ball_lambda = B.line_positive_intersect(r, v);
        int facet = P.num_of_hyperplanes();

        if (polypair.first < ball_lambda ) facet = polypair.second;

        return std::pair<NT, int>(std::min(polypair.first, ball_lambda), facet);
    }


    //First coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        std::pair <NT, NT> polypair = P.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <NT, NT> ballpair = B.line_intersect_coord(r, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    //Not the first coordinate ray shooting intersecting convex body
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        std::pair <NT, NT> polypair = P.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        std::pair <NT, NT> ballpair = B.line_intersect_coord(r, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<NT,NT> query_dual(Point &p, unsigned int rand_coord) {
        std::pair <NT, NT> polypair = P.query_dual(p, rand_coord);
        std::pair <NT, NT> ballpair = B.line_intersect_coord(p, rand_coord);
        return std::pair<NT, NT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    void compute_reflection (Point &v, Point &p, int &facet) {

        if (facet == P.num_of_hyperplanes()) {
            B.compute_reflection(v, p, facet);
        } else {
            P.compute_reflection(v, p, facet);
        }

    }

    void shift(VT et){}

    void normalize() {
        P. normaize();
    }

    template <class T>
    bool get_points_for_rounding (T &randPoints) {
        return false;
    }

    std::vector<NT> get_dists(NT radius){
        return P.get_dists(radius);
    }

    void linear_transformIt(MT T) {}

    void free_them_all(){
        P.free_them_all();
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
