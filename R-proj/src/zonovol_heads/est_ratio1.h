// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef EST_RATIO_ONE_H
#define EST_RATIO_ONE_H

template <class Zonotope, class HPolytope, typename NT, class Parameters>
NT est_ratio_hzono(Zonotope &Z, HPolytope &HP, NT error, NT ratio, Parameters &var, NT &steps) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    bool print = var.verbose;

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    int W=4*n*n+500;
    int m = Z.num_of_generators();
    NT curr_eps = error;
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
    NT countIn = ratio*NT(1200);
    NT totCount = NT(1200);
    Point p(n);
    while(!done){

        //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
        //sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
        rand_point(HP, p, var);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
            //q2=*rpit;
        }
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
            if(print) std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            steps = (totCount - 1200.0);
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;

    }
    std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
    return val;
}



template <class Point, class convexB, class Zonotope, class Parameters, typename NT>
NT est_ratio_zonoballs(Zonotope &Z, convexB &b1, NT ratio, NT error, Parameters &var, NT &steps){


    //typedef typename Zonotope::PolytopePoint Point;
    //ball b1 = zb1.second();
   // ball b2 = zb2.second();

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    bool print = var.verbose;

    //typedef typename Zonotope::PolytopePoint Point;
    int n = var.n;
    NT countIn = ratio*NT(1200+2.0*n*n);
    NT totCount = NT(1200+2.0*n*n);
    int W=4*n*n+500;
    //int m = Z.num_of_generators();
    NT curr_eps = error;
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
    NT val = ratio;

    // Set the radius for the ball walk if it is requested
    //if (var.ball_walk) {
    //if (var.deltaset) {
    //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
    //}
    //}

    //NT countIn = 0.0;
    //NT totCount = 0.0;
    Point p(n);
    while(!done){

        rand_point(Z, p, var);
        //Point q2;
        //for ( ;  rpit!=randPoints.end(); ++rpit) {
        if (b1.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }

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
            if(print) std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            steps = (totCount - 1200.0);
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;
        //}
    }

}

template <typename NT>
bool check_max_error2(NT a, NT b, NT val, NT error) {

    NT e1 = std::abs(a - val) / a;
    NT e2 = std::abs(b - val) / b;
    //std::cout<<"er1 = "<<e1<<" er2 = "<<e2<<" error = "<<error<<std::endl;
    if (e1<error/2.0 && e2<error/2.0){
        return true;
    }
    return false;
}


template <class Zonotope, class HPolytope, typename NT, class Parameters>
NT est_ratio_hzono_normal(Zonotope &Z, HPolytope &HP, NT error, int n_subw, int n_tuple, NT prob, NT ratio, Parameters &var, NT &steps) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    //int W = 4 * n * n + 500;
    //bool print = var.verbose;
    //int m = Z.num_of_generators();
    int W = n_subw*n_tuple;
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
    std::vector<std::vector<NT> > vals(n_subw, std::vector<NT>(n_tuple));

    NT countIn = ratio*1200.0;
    NT totCount = 1200.0;
    //NT rad = B0.radius();
    int col=0, row=0;
    std::vector<NT> sums(n_subw,0.0);
    Point p(n);
    //typename std::vector< std::vector<NT> >::iterator sumit=vals.begin();
    //typename std::vector<NT>::iterator sit=sums.begin();
    for (int i = 0; i < W; ++i) {
        rand_point(HP, p, var);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;
        vals[row][col] = val;
        col++;
        if(col%n_tuple==0) {
            col=0;

            //sums[row-1] = std::accumulate((*sumit).begin(), (*sumit).end(), NT(0.0)) / 120.0;
            for (int j = 0; j < n_tuple; ++j) {
                sums[row] += vals[row][j];
            }
            sums[row] = sums[row]/NT(n_tuple);
            row++;
            //0sumit++;
            //sit++;
        }
    }
    //std::cout<<"countIn = "<<countIn<<" totCount = "<<totCount<<std::endl;
    col=0;
    std::pair<NT,NT> mv;
    NT pr = (1.0 + prob) / 2.0 ,m, s, zp;
    bool chk;
    //p=Point(n);
    while(!done){

        rand_point(HP, p, var);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;

        val = countIn / totCount;
        for (int i = 0; i < n_subw-1; ++i) {
            sums[i] -= vals[i][col] / NT(n_tuple);
            sums[i] += vals[i+1][col] / NT(n_tuple);
            vals[i][col] = vals[i+1][col];
        }
        sums[n_subw-1] -= vals[n_subw-1][col] / NT(n_tuple);
        sums[n_subw-1] += val / NT(n_tuple);
        vals[n_subw-1][col] = val;

        col++;
        if(col%n_tuple==0) col=0;

        mv = getMeanVariance(sums);
        m = mv.first;
        s = std::sqrt(mv.second);
        //zp = m + s*std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
        zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
        chk= check_max_error2(m-zp*s, m+zp*s, val, error);
        if(chk) {
            //std::cout<<"final rejection to Z ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            steps = (totCount - 1200.0);
            return val;
        }

    }
    return val;
}



#endif
