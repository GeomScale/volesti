// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef CG_ZONOVOL_H
#define CG_ZONOVOL_H

template <class Polytope, class UParameters, class GParameters, class Point, typename NT, class VT, class MT>
NT cg_volume_zono(Polytope &P,
                             GParameters &var,  // constans for volume
                             UParameters &var2,
                             std::pair<Point,NT> InnerBall,
                             Rcpp::Function rtmvnorm, Rcpp::Function mvrandn, Rcpp::Function mvNcdf,
                             MT sigma, VT l, VT u, int Wst) {
    //typedef typename Polytope::MT 	MT;
    //typedef typename Polytope::VT 	VT;
    typedef typename UParameters::RNGType RNGType;
    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    NT vol;
    bool round = var.round, done;
    bool print = var.verbose;
    bool rand_only = var.rand_only, deltaset = false;
    unsigned int n = var.n, steps;
    unsigned int walk_len = var.walk_steps, m=P.num_of_hyperplanes();
    unsigned int n_threads = var.n_threads, min_index, max_index, index, min_steps;
    NT error = var.error, curr_eps, min_val, max_val, val;
    NT frac = var.frac;
    RNGType &rng = var.rng;
    typedef typename std::vector<NT>::iterator viterator;

    // Consider Chebychev center as an internal point
    Point c=InnerBall.first;
    NT radius=InnerBall.second;
    if (var.ball_walk){
        if(var.delta<0.0){
            var.delta = 4.0 * radius / NT(n);
            var.deltaset = true;
        }
    }

    // rounding of the polytope if round=true
    NT round_value=1;
    if(round){
#ifdef VOLESTI_DEBUG
        if(print) std::cout<<"\nRounding.."<<std::endl;
#endif
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::pair<NT,NT> res_round = rounding_min_ellipsoid(P,InnerBall,var2);
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
#ifdef VOLESTI_DEBUG
        if(print) std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
#endif
        round_value=res_round.first;
        std::pair<Point,NT> res=P.ComputeInnerBall();
        c=res.first; radius=res.second;
    }

    // Save the radius of the Chebychev ball
    var.che_rad = radius;

    // Move chebychev center to origin and apply the same shifting to the polytope
    //VT c_e(n);
    //for(unsigned int i=0; i<n; i++){
        //c_e(i)=c[i];  // write chebychev center in an eigen vector
    //}
    //P.shift(c_e);
    MT G = P.get_mat().transpose();

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = var.ratio;
    NT C = var.C;
    unsigned int N = var.N;
    std::vector<NT> probs;

    // Computing the sequence of gaussians
#ifdef VOLESTI_DEBUG
    if(print) std::cout<<"\n\nComputing annealing...\n"<<std::endl;
#endif
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
    get_annealing_schedule2(P, radius, ratio, C, frac, N, var, error, a_vals, Wst, l, u, sigma, rtmvnorm, mvrandn, mvNcdf, probs);
    double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
#ifdef VOLESTI_DEBUG
    if(print) std::cout<<"All the variances of schedule_annealing computed in = "<<tstop2-tstart2<<" sec"<<std::endl;
#endif

    unsigned int mm = a_vals.size()-1, j=0;
    if(print){
        for (viterator avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++){
            std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
        }
        std::cout<<"\n"<<std::endl;
    }

    // Initialization for the approximation of the ratios
    std::vector<NT> fn(mm,0), its(mm,0), lamdas(m,0);
    unsigned int W = var.W;
    std::vector<NT> last_W2(W,0);
    vol=std::pow(M_PI/a_vals[0], (NT(n))/2.0)*std::abs(round_value);
    Point p(n); // The origin is in the Chebychev center of the Polytope
    Point p_prev=p;
    unsigned int coord_prev, i=0;
    viterator fnIt = fn.begin(), itsIt = its.begin(), avalsIt = a_vals.begin(), probIt = probs.begin(), minmaxIt;

#ifdef VOLESTI_DEBUG
    if(print) std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    if(print) std::cout<<"computing ratios..\n"<<std::endl;
#endif

    // Compute the first point if CDHR is requested
    //if(var.coordinate && !var.ball_walk){
        //gaussian_first_coord_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
    //}
    MT pointset;
    int kk = P.num_of_generators();
    double tstart122 = (double)clock()/(double)CLOCKS_PER_SEC;
    for ( ; fnIt != fn.end(); fnIt++, itsIt++, avalsIt++, ++probIt, i++) { //iterate over the number of ratios
        //initialize convergence test
        curr_eps = error/std::sqrt((NT(mm)));
        done=false;
        min_val = minNT;
        max_val = maxNT;
        min_index = W-1;
        max_index = W-1;
        index = 0;
        min_steps=0;
        std::vector<NT> last_W=last_W2;

        // Set the radius for the ball walk if it is requested
        //if (var.ball_walk) {
            //if (var.deltaset) {
                //var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
            //}
        //}

                MT sigma2;
        while(!done || (*itsIt)<min_steps){

            //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
            sigma2 = (1.0/(2.0*(*avalsIt)))*sigma;
            if((*probIt)>0.001) {
                std::cout<<"a_"<<i<<" = "<<*avalsIt<<" | prob = "<<*probIt<<std::endl;
                pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);
            } else {
                std::cout<<"a_"<<i<<" = "<<*avalsIt<<" | prob = "<<*probIt<<std::endl;
                pointset = sampleTr_gibbs(l, u, sigma2, 2*W,  Wst, rtmvnorm, G);
            }
            //pointset = sampleTr(l, u , sigma2, 2*W, mvrandn, G);

            for (int k = 0; k < 2*W; ++k) {
                *itsIt = *itsIt + 1.0;
                //*fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
                *fnIt = *fnIt + std::exp(-(*(avalsIt + 1))*(pointset.col(k).squaredNorm())) / std::exp(-(*avalsIt)*(pointset.col(k).squaredNorm()));
                val = (*fnIt) / (*itsIt);

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
                    done=true;
                }

                index = index%W+1;

                if(index==W) index=0;
            }
        }
#ifdef VOLESTI_DEBUG
        if(print) std::cout<<"ratio "<<i<<" = "<<(*fnIt) / (*itsIt)<<" N_"<<i<<" = "<<*itsIt<<std::endl;
#endif
        vol = vol*((*fnIt) / (*itsIt));
    }
    double tstop122 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(print) std::cout<<"Ratios time cost = "<<tstop122-tstart122<<" sec"<<std::endl;
    // Compute and print total number of steps in verbose mode only

    //if (print) {
    NT sum_of_steps = 0.0;
    for(viterator it = its.begin(); it != its.end(); ++it) {
        sum_of_steps += *it;
    }
    //cg_steps = sum_of_steps + (NT(N))*(NT(mm));
    steps= int(sum_of_steps);
#ifdef VOLESTI_DEBUG
       if(print) std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
#endif
    //}


    return vol;
}

#endif
