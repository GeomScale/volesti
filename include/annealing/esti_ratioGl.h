// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef ESTI_RATIOGL_H
#define ESTI_RATIOGL_H

#include <numeric>

template <class ZonoBall, class ball, typename NT, class Parameters>
NT esti_ratio(ZonoBall &Zb, ball B0, NT ratio, NT error, int Win, Parameters &var, NT &steps) {

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
    if(var.coordinate && !var.ball_walk){
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
            steps = (totCount - 1200.0-n*n*2.0);
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
NT esti_ratio2(ball B0, Zonotope &Z, NT error, int Win, NT ratio, NT &steps) {

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
            steps = (totCount - 1200.0);
            return val;
        }

        index = index%W+1;

        if(index==W) index=0;

    }
    return countIn / totCount;
}

template <typename NT>
bool check_max_error123(NT a, NT b, NT val, NT error) {

    NT e1 = std::abs(a - val) / a;
    NT e2 = std::abs(b - val) / b;
    NT e3 = (b-a)/b;
    //std::cout<<"er1 = "<<e1<<" er2 = "<<e2<<"e3 = "<<e3<<" error = "<<error/10.0<<std::endl;
    //if (e1<error/2.0 && e2<error/2.0){
    if(e3<error/2.0) {
        return true;
    }
    return false;
}

template <class Point, class ZonoBall, class ball, typename NT, class Parameters>
NT esti_ratio_interval(ZonoBall &Zb, ball B0, NT ratio, NT error, int WW, NT prob, Parameters &var, NT &steps) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    //typedef typename ball::BallPoint Point;
    int n = var.n;
    bool print = var.verbose;
    //std::cout<<"n = "<<n<<std::endl;
    //int W=4*n*n+500;
    int W = WW;
    //int m = Z.num_of_generators();
    //std::cout<<"W = "<<W<<" walk_steps = "<<var.walk_steps<<std::endl;
    //std::cout<<"coordinate : "<<var.coordinate<<std::endl;
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
    //std::vector<std::vector<NT> > vals(n_subw, std::vector<NT>(n_tuples));

    NT countIn = ratio*(1200.0+2.0*n*n);
    NT totCount = 1200.0+n*n*2.0;
    //if (print) std::cout<<"countIn = "<<countIn<<" totCount = "<<totCount<<std::endl;
    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;
    if(var.coordinate && !var.ball_walk){
        uniform_first_coord_point(Zb,p,p_prev,coord_prev,var.walk_steps,lamdas,var);
    }
    int col=0, row=0;
    //std::vector<NT> sums(n_subw,0.0);
    //typename std::vector< std::vector<NT> >::iterator sumit=vals.begin();
    //typename std::vector<NT>::iterator sit=sums.begin();
    NT sum_sq=0.0;
    NT sum=0.0;
    for (int i = 0; i < W; ++i) {
        //std::cout<<"row = "<<row<<" col = "<<col<<"i = "<<i<<std::endl;
        uniform_next_point(Zb, p, p_prev, coord_prev, 1, lamdas, var);
        if(B0.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;
        sum += val;
        sum_sq += val*val;
        last_W[index] = val;
        index = index%W+1;

        if(index==W) index=0;
        /*vals[row][col] = val;
        col++;
        if(col%n_tuples==0) {
            col=0;

            //sums[row-1] = std::accumulate((*sumit).begin(), (*sumit).end(), NT(0.0)) / 120.0;
            for (int j = 0; j < n_tuples; ++j) {
                sums[row] += vals[row][j];
            }
            sums[row] = sums[row]/NT(n_tuples);
            row++;
            //0sumit++;
            //sit++;
        }*/
    }
    col=0;
    std::pair<NT,NT> mv;
    NT pr = (1.0 + prob) / 2.0 ,m, mW, s, zp, zp2=2.38774;
    zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
    //std::cout<<"zp = "<<zp<<std::endl;
    m=sum/NT(W);
    NT m2,s2;
    bool chk;
    while(!done){
        //std::cout<<"col = "<<col<<std::endl;

        uniform_next_point(Zb, p, p_prev, coord_prev, 1, lamdas, var);
        if(B0.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;


        /*for (int i = 0; i < n_subw-1; ++i) {
            sums[i] -= vals[i][col] / NT(n_tuples);
            sums[i] += vals[i+1][col] / NT(n_tuples);
            vals[i][col] = vals[i+1][col];
        }
        sums[n_subw-1] -= vals[n_subw-1][col] / NT(n_tuples);
        sums[n_subw-1] += val / NT(n_tuples);
        vals[n_subw-1][col] = val;

        col++;
        if(col%n_tuples==0) col=0;

        //mv = getMeanVariance(sums);
        //m = mv.first;
        //s = std::sqrt(mv.second);*/
        m-=last_W[index]/NT(W);
        m+=val/NT(W);
        sum_sq -= last_W[index]*last_W[index];
        sum_sq += val*val;
        sum-=last_W[index];
        sum+=val;
        //mv = getMeanVariance(last_W);
        //m2=mv.first;
        s = std::sqrt((sum_sq+NT(W)*m*m-2.0*m*sum)/NT(W));
        //mv = getMeanVariance(last_W);
        //m2=mv.first;
        //s2 = std::sqrt(mv.second);
        //std::cout<<"m="<<m<<" m2="<<m2<<"s="<<s<<" s2="<<s2<<std::endl;
        last_W[index] = val;
        index = index%W+1;

        if(index==W) index=0;
        //zp2 = m + s*std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);

        //std::cout<<"["<<m-zp*s<<","<<m+zp*s<<"] || ["<<m-zp2*s<<","<<m+zp2*s<<"]"<<std::endl;
        chk= check_max_error123(m-zp*s, m+zp*s, val, error);
        if(chk) {
            if (print) std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
            done=true;
            steps = (totCount - 1200.0-n*n*2.0)*NT(var.walk_steps);
            //std::cout<<"steps = "<<steps<<std::endl;
            return val;
        }

    }
    return val;

}

template <class RNGType, class Zonotope, class ball, typename NT>
NT esti_ratio2_const(ball B0, Zonotope &Z, NT error, int WW, NT ratio, NT prob, NT &steps) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    //int W = 4 * n * n + 500;
    //bool print = var.verbose;
    //int m = Z.num_of_generators();
    int W = WW;
   // std::cout<<"W = "<<W<<" n_subw = "<<n_subw<<" n_tuples = "<<n_tuple<<std::endl;
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
    //std::vector<std::vector<NT> > vals(n_subw, std::vector<NT>(n_tuple));

    NT countIn = ratio*1200.0;
    NT totCount = 1200.0;
    NT rad = B0.radius();
    int col=0, row=0;
    //std::vector<NT> sums(n_subw,0.0);
    Point p(n);
    //typename std::vector< std::vector<NT> >::iterator sumit=vals.begin();
    //typename std::vector<NT>::iterator sit=sums.begin();
    NT sum_sq=0.0;
    NT sum=0.0;
    for (int i = 0; i < W; ++i) {
        p = get_point_in_Dsphere<RNGType, Point>(n, rad);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
        val = countIn / totCount;
        sum += val;
        sum_sq += val*val;
        last_W[index] = val;
        index = index%W+1;

        if(index==W) index=0;
        /*vals[row][col] = val;
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
        }*/
    }
    //std::cout<<"countIn = "<<countIn<<" totCount = "<<totCount<<std::endl;
    col=0;
    std::pair<NT,NT> mv;
    NT pr = (1.0 + prob) / 2.0 ,m, s, zp;
    zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);
   // std::cout<<"zp = "<<zp<<std::endl;
    bool chk;
    //mv = getMeanVariance(last_W);
    m=sum/NT(W);
    //s = std::sqrt(mv.second);

    while(!done){

        p = get_point_in_Dsphere<RNGType, Point>(n, rad);
        if(Z.is_in(p)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;

        val = countIn / totCount;

        /*for (int i = 0; i < n_subw-1; ++i) {
            sums[i] -= vals[i][col] / NT(n_tuple);
            sums[i] += vals[i+1][col] / NT(n_tuple);
            vals[i][col] = vals[i+1][col];
        }
        sums[n_subw-1] -= vals[n_subw-1][col] / NT(n_tuple);
        sums[n_subw-1] += val / NT(n_tuple);
        vals[n_subw-1][col] = val;

        col++;
        if(col%n_tuple==0) col=0;

        //mv = getMeanVariance(sums);
        //m = mv.first;
        //s = std::sqrt(mv.second);*/
        m-=last_W[index]/NT(W);
        m+=val/NT(W);
        sum_sq -= last_W[index]*last_W[index];
        sum_sq += val*val;
        sum-=last_W[index];
        sum+=val;
        //mv = getMeanVariance(last_W);
        //m=mv.first;
        s = std::sqrt((sum_sq+NT(W)*m*m-2.0*m*sum)/NT(W));

        last_W[index] = val;
        index = index%W+1;

        if(index==W) index=0;
        //zp = m + s*std::sqrt(2.0)*boost::math::erf_inv(2.0*pr - 1.0);

        chk= check_max_error123(m-zp*s, m+zp*s, val, error);
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
