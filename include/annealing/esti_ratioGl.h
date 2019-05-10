// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef ESTI_RATIOGL_H
#define ESTI_RATIOGL_H

#include <numeric>
//#include <boost/math/distributions/students_t.hpp>

template <class ZonoBall, class ball, typename NT, class Parameters>
NT esti_ratio(ZonoBall &Zb, ball B0, NT ratio, NT error, int Win, Parameters &var) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    typedef typename ball::BallPoint Point;
    int n = var.n;
    bool print = var.verbose;
    //std::cout<<"n = "<<n<<std::endl;
    int W=4*n*n+500;
    //int W = Win;
    //int m = Z.num_of_generators();
    NT curr_eps = error;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0), lamdas(Zb.num_of_hyperplanes(),0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    NT val;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    //MT sigma2;
    //MT sample;
    NT countIn = ratio*(1200.0+2.0*n*n);
    NT totCount = 1200.0+n*n*2.0;
    //if (print) std::cout<<"countIn = "<<countIn<<" totCount = "<<totCount<<std::endl;
    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;
    if(var.cdhr_walk && !var.ball_walk){
        uniform_first_coord_point(Zb,p,p_prev,coord_prev,var.walk_steps,lamdas,var);
    }
    while(!done){

        //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
        //sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
        //rand_point(Zb, p, var);
        uniform_next_point(Zb, p, p_prev, coord_prev, var.walk_steps, lamdas, var);
        if(B0.is_in(p)==-1) {
            //std::cout<<"Point in!"<<std::endl;
            countIn = countIn + 1.0;
            //q2=*rpit;
        }
        //std::cout<<"Point out!"<<std::endl;
        totCount = totCount + 1.0;

        //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

        //for (int k = 0; k < 4*W; ++k) {
        //*itsIt = *itsIt + 1.0;
        //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
        // *fnIt = *fnIt + std::exp(-(*(avalsIt + 1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
        //val = (*fnIt) / (*itsIt);

        val = countIn / totCount;
        last_W[index] = val;
        if(val<=min_val){
            min_val = val;
            min_index = index;
        }else if(min_index==index){
            minmaxIt = std::min_element(last_W.begin(), last_W.end());
            min_val = *minmaxIt;
            min_index = std::distance(last_W.begin(), minmaxIt);
        }

        if(val>=max_val){
            max_val = val;
            max_index = index;
        }else if(max_index==index){
            minmaxIt = std::max_element(last_W.begin(), last_W.end());
            max_val = *minmaxIt;
            max_index = std::distance(last_W.begin(), minmaxIt);
        }

        if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
            if (print) std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            //steps = (totCount - 1200.0-n*n*2.0);
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;

    }
    std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
    return val;
}


template <class RNGType, class Zonotope, class ball, typename NT>
NT esti_ratio1(ball B0, Zonotope &Z, NT ratio, int N) {

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    int W = 10 * n * n + 1200;
    int m = Z.num_of_generators();

    NT countIn = ratio * 1200.0;
    NT totCount = 1200.0;
    NT rad = B0.radius();
    Point p(n);
    for (int i = 0; i < N; ++i) {
        p = get_point_in_Dsphere<RNGType, Point>(n, rad);
        if (Z.is_in(p) == -1) {
            countIn = countIn + 1.0;
            //q2=*rpit;
        }
        totCount = totCount + 1.0;
    }
    return countIn / totCount;
}

template <class RNGType, class Zonotope, class ball, typename NT>
NT esti_ratio2(ball B0, Zonotope &Z, NT error, int Win, NT ratio) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    int W = 4 * n * n + 500;
    //bool print = var.verbose;
    //int m = Z.num_of_generators();
    //int W = Win;
    // std::cout<<"W = "<<W<<std::endl;
    NT curr_eps = error;
    //std::cout<<"curr_eps = "<<curr_eps<<std::endl;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    NT val;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    //MT sigma2;
    //MT sample;
    NT countIn = ratio*1200.0;
    NT totCount = 1200.0;
    NT rad = B0.radius();
    //std::cout<<"countIn = "<<countIn<<" totCount = "<<totCount<<std::endl;
    Point p(n);
    while(!done){

        //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
        //sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
        p = get_point_in_Dsphere<RNGType, Point>(n, rad);
        if(Z.is_in(p)==-1) {
            //std::cout<<"Point in!"<<std::endl;
            countIn = countIn + 1.0;
            //q2=*rpit;
        }
        //std::cout<<"Point out!"<<std::endl;
        totCount = totCount + 1.0;

        //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

        //for (int k = 0; k < 4*W; ++k) {
        //*itsIt = *itsIt + 1.0;
        //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
        // *fnIt = *fnIt + std::exp(-(*(avalsIt +  1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
        //val = (*fnIt) / (*itsIt);

        val = countIn / totCount;
        last_W[index] = val;
        if(val<=min_val){
            min_val = val;
            min_index = index;
        }else if(min_index==index){
            minmaxIt = std::min_element(last_W.begin(), last_W.end());
            min_val = *minmaxIt;
            min_index = std::distance(last_W.begin(), minmaxIt);
        }

        if(val>=max_val){
            max_val = val;
            max_index = index;
        }else if(max_index==index){
            minmaxIt = std::max_element(last_W.begin(), last_W.end());
            max_val = *minmaxIt;
            max_index = std::distance(last_W.begin(), minmaxIt);
        }

        if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
            //std::cout<<"last ball rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
           // steps = (totCount - 1200.0);
            //std::cout<<"COUNT IN Cm = "<<steps<<std::endl;
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;

    }
    return countIn / totCount;
}

template <typename NT>
bool check_max_error123(NT a, NT b, NT val, NT error) {

    //NT e1 = std::abs(a - val) / a;
    //NT e2 = std::abs(b - val) / b;
    NT e3 = (b-a)/a;
    if(e3<error/2.0) {
        return true;
    }
    return false;
}

template <class RNGType, class Point, class PolyBall1, class PolyBall2, typename NT, class Parameters>
NT esti_ratio_interval(PolyBall1 &Pb1, PolyBall2 Pb2, NT ratio, NT error, int W, int Ntot, NT prob,
                       Parameters &var, bool isball = false, NT radius = 0.0) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    int n = var.n;
    bool print = var.verbose;
    NT curr_eps = error;
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    std::vector<NT> last_W(W,0), lamdas(Pb1.num_of_hyperplanes(),0);
    NT val;

    NT countIn = ratio*NT(Ntot);
    NT totCount = NT(Ntot);
    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;
    if(var.cdhr_walk && !var.ball_walk && !isball){
        uniform_first_coord_point(Pb1,p,p_prev,coord_prev,var.walk_steps,lamdas,var);
    }
    int col=0, row=0;
    NT sum_sq=0.0;
    NT sum=0.0;
    for (int i = 0; i < W; ++i) {

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, 1, lamdas, var);
        }
        if (Pb2.is_in(p) == -1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;
        sum += val;
        sum_sq += val * val;
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
    }

    NT pr = (1.0 + prob) / 2.0 ,m, mW, s, zp, zp2=2.38774;
    zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
    m=sum/NT(W);

    while(!done) {

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, 1, lamdas, var);
        }
        if (Pb2.is_in(p) == -1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;

        m -= last_W[index] / NT(W);
        m += val / NT(W);
        sum_sq -= last_W[index] * last_W[index];
        sum_sq += val * val;
        sum -= last_W[index];
        sum += val;
        s = std::sqrt((sum_sq + NT(W) * m * m - 2.0 * m * sum) / NT(W));
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
        if (check_max_error123(val - zp * s, val + zp * s, val, error)) {
            if (print) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
            //done=true;
            return val;
        }

    }
    return val;

}

/*
template <class RNGType, class Zonotope, class ball, typename NT>
NT esti_ratio2_const(ball B0, Zonotope &Z, NT error, int W, NT ratio, NT prob) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    bool done=false;
    NT min_val = minNT;
    NT max_val = maxNT;
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int min_steps=0;
    std::vector<NT> last_W(W,0);
    NT val;

    NT countIn = ratio*1200.0;
    NT totCount = 1200.0;
    NT rad = B0.radius();
    int col=0, row=0;
    Point p(n);
    NT sum_sq=0.0;
    NT sum=0.0;
    for (int i = 0; i < W; ++i) {
        p = get_point_in_Dsphere<RNGType, Point>(n, rad);
        if (Z.is_in(p) == -1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;
        sum += val;
        sum_sq += val * val;
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
    }
    std::pair<NT,NT> mv;
    NT pr = (1.0 + prob) / 2.0 ,m, s, zp;
    zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
    bool chk;
    m=sum/NT(W);

    while(!done){

        p = get_point_in_Dsphere<RNGType, Point>(n, rad);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;

        val = countIn / totCount;
        m-=last_W[index]/NT(W);
        m+=val/NT(W);
        sum_sq -= last_W[index]*last_W[index];
        sum_sq += val*val;
        sum-=last_W[index];
        sum+=val;
        s = std::sqrt((sum_sq+NT(W)*m*m-2.0*m*sum)/NT(W));

        last_W[index] = val;
        index = index%W+1;

        if(index==W) index=0;
        chk= check_max_error123(val-zp*s, val+zp*s, val, error);
        if(chk) {
            std::cout<<"final rejection to Z ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            //done=true;
            return val;
        }

    }
    return val;
}*/

/*
template <class RNGType, class Zonotope, class ball, typename NT>
NT esti_ball_ratio_no_W(ball B0, Zonotope &Z, NT error, std::vector<NT> all_ratios,
                     NT ratio, NT prob, NT &steps, int N) {

    N=120;
    NT mean = ratio, s, zp, M2 = 0, variance = 0, delta, pr = (1.0 + prob) / 2.0, rad = B0.radius(), count_in,
            tot_in = NT(N), alpha = 1.0-prob, T;
    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension(), nu = 0;
    Point p(n);
    zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
    //std::pair<NT,NT> mv = getMeanVariance(all_ratios);

    typedef typename std::vector<NT>::iterator viterator;
    viterator vecit = all_ratios.begin();
    for( ; vecit!=all_ratios.end(); vecit++, nu++){
        std::cout<<*vecit<<" "<<std::endl;
        delta = *vecit - mean;
        mean += delta / (nu + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (nu + 1);
    }
    std::cout<<"\nmean = "<<mean<<" ratio = "<<ratio<<" std = "<<std::sqrt(variance)
             <<" nu = "<<nu<<" N = "<<N<<" zp = "<<zp<<" error = "<<error
             <<" er ="<<(zp*std::sqrt(variance))/(mean - zp*std::sqrt(variance))<<std::endl;

    int count = 0;
    while(true){

        count_in = 0.0;
        for (int i = 0; i < N; ++i) {
            count++;
            p = get_point_in_Dsphere<RNGType, Point>(n, rad);
            if (Z.is_in(p)==-1) {
                count_in = count_in + 1.0;
            }
        }
        delta = count_in/tot_in - mean;
        mean += delta / (nu + 1);
        M2 += delta * (count_in/tot_in - mean);
        variance = M2 / (nu + 1);
        nu++;
        s = std::sqrt(variance);

        boost::math::students_t dist(nu - 1);
        T = boost::math::quantile(boost::math::complement(dist, alpha / 2.0))/std::sqrt(NT(nu));

        //std::cout<<"ratio = "<<count_in/tot_in<<" mean = "<<mean<<" std = "<<s
               //  <<" er = "<<(T*s)/(mean - T*s)<<std::endl;

        if ((T*s)/(mean - T*s) <= error) {
            std::cout<<"COUNT IN Cm = "<<count<<" nu = "<<nu<<std::endl;
            return mean;
        }

    }

}*/

#endif
