// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef EST_RATIO_ONE_H
#define EST_RATIO_ONE_H

template <class Zonotope, class MT, class VT, typename NT>
NT est_ratio_zono(Zonotope &Z, NT &prob, NT error, int W, Rcpp::Function rtmvnorm,
                    Rcpp::Function mvrandn, MT &sigma, MT &G, VT &l, VT &u) {

    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    W=10*n*n+1200;
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

    MT sigma2;
    MT sample;
    NT countIn = 0.0;
    NT totCount = 0.0;
    while(!done){

        //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
        //sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
        if(prob>0.001) {
            sample = sampleTr(l, u , sigma, 4*W, mvrandn, G);
        } else {
            sample = sampleTr_gibbs(l, u, sigma, 4*W, 20+m/5, rtmvnorm, G);
        }
        randPoints.clear();
        for (int i = 0; i < 4*W; ++i) {
            Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
            randPoints.push_back(p);
        }
        rpit = randPoints.begin();

        //Point q2;
        for ( ;  rpit!=randPoints.end(); ++rpit) {
            if(Z.is_in(*rpit)==-1) {
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
                std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
                done=true;
                return val;
            }

            index = index%W+1;

            if(index==W) index=0;
        }
    }
    std::cout<<"final rejection ratio = "<<val<< " | total points = "<<totCount<<std::endl;
    return val;
}

#endif
