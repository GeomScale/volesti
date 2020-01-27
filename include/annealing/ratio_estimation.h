// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef RATIO_ESTIMATION_H
#define RATIO_ESTIMATION_H

#define MAX_ITER_ESTI 500000

template <typename NT>
bool check_max_error(const NT &a, const NT &b, const NT &error) {

    if((b-a)/a<error/2.0) {
        return true;
    }
    return false;

}


template <typename RNGType, typename Point, typename PolyBall1, typename PolyBall2, typename NT, typename Parameters>
NT esti_ratio(PolyBall1 &Pb1, PolyBall2 &Pb2, const NT &ratio, const NT &error, const int &W,
        const int &Ntot, const Parameters &var, bool isball = false, NT radius = 0.0) {

    int n = var.n, min_index = W-1, max_index = W-1, index = 0, iter = 1;
    bool print = var.verbose;
    NT min_val = std::numeric_limits<NT>::lowest(), max_val = std::numeric_limits<NT>::max(), val, lambda;
    size_t totCount = Ntot, countIn = Ntot * ratio;
    std::vector<NT> last_W(W), lamdas(Pb1.num_of_hyperplanes()), Av(Pb1.num_of_hyperplanes());
    std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    typename std::list<Point>::iterator rpit;
    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;

    if(!var.ball_walk && !isball){
        uniform_first_point(Pb1,p,p_prev,coord_prev,var.walk_steps,lamdas,Av,lambda,var);
    }

    while(iter <= MAX_ITER_ESTI){
        iter++;

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps, lamdas, Av, lambda, var);
        }
        if(Pb2.is_in(p)==-1) countIn = countIn + 1.0;

        totCount = totCount + 1.0;
        val = NT(countIn) / NT(totCount);
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
            return val;
        }

        index = index%W+1;
        if(index==W) index=0;

    }
    return val;
}


template <typename RNGType, typename Point, typename PolyBall1, typename PolyBall2, typename NT, typename Parameters>
NT esti_ratio_interval(PolyBall1 &Pb1, PolyBall2 &Pb2, const NT &ratio, const NT &error, const int &W,
        const int &Ntot, const NT &prob, const Parameters &var, bool isball = false, NT radius = 0.0) {

    int n = var.n, index = 0, iter = 1;
    bool print = var.verbose;
    std::vector<NT> last_W(W), lamdas(Pb1.num_of_hyperplanes()), Av(Pb1.num_of_hyperplanes());
    NT val, sum_sq=0.0, sum=0.0, lambda;
    size_t totCount = Ntot, countIn = Ntot * ratio;
    //std::cout<<"countIn = "<<countIn<<", totCount = "<<totCount<<std::endl;

    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;
    if(!var.ball_walk && !isball) uniform_first_point(Pb1, p, p_prev, coord_prev, 1,
                                                                             lamdas, Av, lambda, var);
    for (int i = 0; i < W; ++i) {

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps, lamdas, Av, lambda, var);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1;

        totCount = totCount + 1;
        val = NT(countIn) / NT(totCount);
        sum += val;
        sum_sq += val * val;
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
    }

    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0)), m=sum/NT(W), s;

    while(iter <= MAX_ITER_ESTI) {
        iter++;

        if (isball) {
            p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        } else {
            uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps, lamdas, Av, lambda, var);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1;

        totCount = totCount + 1;
        val = NT(countIn) / NT(totCount);

        m = (m - last_W[index] / NT(W)) + val / NT(W);
        sum_sq = (sum_sq - last_W[index] * last_W[index]) + val * val;
        sum = (sum - last_W[index]) + val;
        s = std::sqrt((sum_sq + NT(W) * m * m - 2.0 * m * sum) / NT(W));
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
        if (check_max_error(val - zp * s, val + zp * s, error)) {
            //if (print) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
            return val;
        }

    }
    return val;

}

#endif

