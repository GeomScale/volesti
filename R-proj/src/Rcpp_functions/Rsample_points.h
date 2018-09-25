// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

#ifndef RSAMPLE_POINTS_H
#define RSAMPLE_POINTS_H

Rcpp::NumericMatrix Rsample_points (Rcpp::NumericMatrix A, unsigned int walk_len, double e, Rcpp::NumericVector InnerVec,
                           bool CG, bool ball_walk, double delta, bool Vpoly, bool Zono, bool sam_simplex,
                           bool sam_can_simplex, bool sam_arb_simplex, bool sam_ball, bool sam_sphere,
                           unsigned int numpoints, double variance) {

    std::list<Point> randPoints;
    Rcpp::NumericMatrix PointSet(n,numpoints);

    if (sam_ball || sam_sphere) {

        for (unsigned int k = 0; k < numpoints; ++k) {
            if (sam_ball) {
                randPoints.push_back(get_point_in_Dsphere<RNGType , Point >(dim_gen, delta));
            } else {
                randPoints.push_back(get_point_on_Dsphere<RNGType , Point >(dim_gen, delta));
            }
        }

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        j = 0;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;

    }

    if (sam_simplex || sam_can_simplex) {

        if (sam_simplex) {
            Sam_Unit<NT, RNGType >(dim_gen, numpoints, randPoints);
        } else {
            Sam_Canon_Unit<NT, RNGType >(dim_gen, numpoints, randPoints);
        }

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        j = 0;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    if (sam_arb_simplex) {
        std::vector<NT> temp_p(n + 1, 0.0);
        typename std::vector<NT>::iterator temp_it;
        std::vector<Point> vec_point;

        for (int k = 0; k < A.nrow(); ++k) {
            temp_it = temp_p.begin();
            for (int l = 0; l < A.ncol(); ++l, ++temp_it) {
                *temp_it = A(k,l);
            }
            vec_point.push_back(Point(n+1, temp_p.begin(), temp_p.end()));
        }

        Sam_arb_simplex<NT, RNGType>(vec_point.begin(), vec_point.end(), numpoints, randPoints);

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        j = 0;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose=false,
            coordinate=coord;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    unsigned int m=A.nrow()-1;
    unsigned int n=A.ncol()-1;
    std::vector <std::vector<NT>> Pin(m + 1, std::vector<NT>(n + 1));

    for (unsigned int i = 0; i < m + 1; i++) {
        for (unsigned int j = 0; j < n + 1; j++) {
            Pin[i][j] = A(i, j);
        }
    }
    // construct polytope
    if (Zono) {
        ZP.init(Pin);
    } else if (!Vpoly) {
        HP.init(Pin);
    } else {
        VP.init(Pin);
    }

    std::pair<Point,NT> InnerBall;
    if (Zono) {
        InnerBall = ZP.ComputeInnerBall();
    } else if (!Vpoly) {
        InnerBall = HP.ComputeInnerBall();
    } else {
        InnerBall = VP.ComputeInnerBall();
    }

    if (InnerVec.size()==n) {
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerVec[j]);
        }
        InnerBall.first = Point( n , temp_p.begin() , temp_p.end() );
    }

    std::list<Point> randPoints;
    Point p = InnerBall.first;
    NT a = 1.0 / (2.0 * variance);
    if (ball_walk) {
        if (delta < 0.0) { // set the radius for the ball walk if is not set by the user
            if (CG) {
                delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
            } else {
                delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
            }
        }
    }
    // initialization
    vars<NT, RNGType> var1(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                           delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
    vars_g<NT, RNGType> var2(n, walk_len, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                             rand_only, false, NN, birk, ball_walk, coord);
    if (Zono) {
        sampling_only<Point>(randPoints, ZP, walk_len, numpoints, CG, a, p, var1, var2);
    } else if (!Vpoly) {
        sampling_only<Point>(randPoints, HP, walk_len, numpoints, CG, a, p, var1, var2);
    } else {
        sampling_only<Point>(randPoints, VP, walk_len, numpoints, CG, a, p, var1, var2);
    }

    // store the sampled points to the output matrix
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    j = 0;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet(i,j)=*qit;
        }
    }
    return PointSet;

}

#endif
