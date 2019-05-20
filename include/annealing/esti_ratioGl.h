// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef ESTI_RATIOGL_H
#define ESTI_RATIOGL_H


template <typename NT>
bool check_max_error123(NT a, NT b, NT error) {

    if((b-a)/a<error/2.0) {
        return true;
    }
    return false;

}


template <class RNGType, class Point, class PolyBall1, class PolyBall2, typename NT, class Parameters>
NT esti_ratio(PolyBall1 &Pb1, PolyBall2 Pb2, NT ratio, NT error, int W, int Ntot, Parameters &var,
              bool isball = false, NT radius = 0.0) {

    int n = var.n, min_index = W-1, max_index = W-1, index = 0;
    bool print = var.verbose;
    NT min_val = std::numeric_limits<NT>::lowest(), max_val = std::numeric_limits<NT>::max(), val,
            countIn = ratio*NT(Ntot), totCount = NT(Ntot);
    std::vector<NT> last_W(W,0), lamdas(Pb1.num_of_hyperplanes(),0);
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;

    if(var.cdhr_walk && !isball){
        uniform_first_coord_point(Pb1,p,p_prev,coord_prev,var.walk_steps,lamdas,var);
    }

    while(true){

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps, lamdas, var);
        }
        if(Pb2.is_in(p)==-1) countIn = countIn + 1.0;

        totCount = totCount + 1.0;
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

        if( (max_val-min_val)/max_val<=error/2.0 ){
            if (print) std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            return val;
        }

        index = index%W+1;
        if(index==W) index=0;
    }
    return val;
}


template <class RNGType, class Point, class PolyBall1, class PolyBall2, typename NT, class Parameters>
NT esti_ratio_interval(PolyBall1 &Pb1, PolyBall2 Pb2, NT ratio, NT error, int W, int Ntot, NT prob,
                       Parameters &var, bool isball = false, NT radius = 0.0) {

    int n = var.n, index = 0;
    bool print = var.verbose;
    std::vector<NT> last_W(W,0), lamdas(Pb1.num_of_hyperplanes(),0);
    NT val, countIn = ratio*NT(Ntot), totCount = NT(Ntot), sum_sq=0.0, sum=0.0;
    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;
    if(var.cdhr_walk && !var.ball_walk && !isball) uniform_first_coord_point(Pb1, p, p_prev, coord_prev, var.walk_steps,
                                                                             lamdas, var);

    for (int i = 0; i < W; ++i) {

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, 1, lamdas, var);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1.0;

        totCount = totCount + 1.0;
        val = countIn / totCount;
        sum += val;
        sum_sq += val * val;
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
    }

    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0)), m=sum/NT(W), s;

    while(true) {

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, 1, lamdas, var);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1.0;

        totCount = totCount + 1.0;
        val = countIn / totCount;

        m = (m - last_W[index] / NT(W)) + val / NT(W);
        sum_sq = (sum_sq - last_W[index] * last_W[index]) + val * val;
        sum = (sum - last_W[index]) + val;
        s = std::sqrt((sum_sq + NT(W) * m * m - 2.0 * m * sum) / NT(W));
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
        if (check_max_error123(val - zp * s, val + zp * s, error)) {
            if (print) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
            return val;
        }

    }
    return val;

}

#endif
