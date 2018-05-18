#include <Rcpp.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <limits>
#include <point_d.h>
#include <polytopes.h>


typedef double                      NT;
typedef point_d								Point;



struct vars{
public:
    vars( int m,
          int n,
          int walk_steps,
          //int n_threads,
          //const double err,
          //const double err_opt,
          //const int lw,
          //double up,
          //const int L,
          //RNGType &rng,
          //generator
          //&get_snd_rand,
          //boost::random::uniform_real_distribution<> urdist,
          //boost::random::uniform_real_distribution<> urdist1,
          bool verbose,
          bool rand_only,
          bool round,
          bool NN,
          bool birk,
          bool coordinate
          ) :
        m(m), n(n), walk_steps(walk_steps), verbose(verbose), rand_only(rand_only), round(round),
        NN(NN),birk(birk),coordinate(coordinate){};

    int m;
    int n;
    int walk_steps;
    //int n_threads;
    //const double err;
   // const double err_opt;
    //const int lw;
    //double up;
    //const int L;
    //RNGType &rng;
   // generator
   // &get_snd_rand;
   // boost::random::uniform_real_distribution<> urdist;
   // boost::random::uniform_real_distribution<> urdist1;
    bool verbose;
    bool rand_only;
    bool round;
    bool NN;
    bool birk;
    bool coordinate;
};




template <class T>
NT volume1_reuse2(T &P,
                  vars &var,  // constans for volume
                  vars &var2, // constants for optimization in case of MinkSums
                  double &Chebtime)
{
    typedef BallIntersectPolytope<T>        BallPoly;

    bool round = var.round;
    bool print = var.verbose;
    bool rand_only = var.rand_only;
    int n = var.n;
    int rnum = var.m;
    int walk_len = var.walk_steps;
    int n_threads = var.n_threads;
    const double err = var.err;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    boost::random::uniform_int_distribution<> uidist(0,n-1);
    //boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

    // Rotation: only for test with skinny polytopes and rounding
    //std::cout<<"Rotate="<<rotate(P)<<std::endl;
    //rotate(P);

    //0. Rounding of the polytope if round=true
    double round_value=1;
    if(round){
        round_value = rounding(P,var,var2);
    }

    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    //1. Compute the Chebychev ball (largest inscribed ball) with center and radius
    double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    if(print) std::cout<<"\nComputing the Chebychev center..."<<std::endl;
    Point c;       //center
    double radius;
    P.chebyshev_center(c,radius);
    //HACK FOR CROSS POLYTOPES
    //std::vector<double> cp(n,0);
    //Point c(n,cp.begin(),cp.end());
    //double radius=std::sqrt(1.0/double(n));

    if(print) std::cout<<"Chebychev center= "<<c<<"\nradius="<<radius<<std::endl;
    double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
    Chebtime = tstop - tstart;
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(print) std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;

    rnum=rnum/n_threads;
    NT vol=0;

    // Perform the procedure for a number of threads and then take the average
    //#pragma omp for ordered schedule(dynamic)
    for(int t=0; t<n_threads; t++){
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        if(print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
        CGAL::Random_points_in_ball_d<Point> gen (n, radius);
        Point p = *gen;
        p = p + (c-CGAL::Origin());
        std::list<Point> randPoints; //ds for storing rand points
        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 50*n, randPoints, var);
        if (print) std::cout<<"First random point: "<<p<<std::endl;

        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        // 3. Sample "rnum" points from P
        if(print) std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
        //randPoints.push_front(p);
        //rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);
        rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        if(print) std::cout << "First random points construction time = " << tstop2 - tstart2 << std::endl;
        //if(rand_only) return -1;

        // 4.  Construct the sequence of balls
        // 4a. compute the radius of the largest ball
        double current_dist, max_dist=NT(0);
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
            current_dist=(*pit-c).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
        if(print) std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist<<std::endl;

        //
        // 4b. Number of balls
        int nb1 = n * (std::log(radius)/std::log(2.0));
        int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));
        //int nb1 = n * (std::log(radius)/std::log(2.0));
        //int nb2 = n * (std::log(max_dist)/std::log(2.0));
        //std::cout<<n* std::log(radius)/std::log(2.0) <<std::endl;
        //std::cout<<n* std::log(max_dist)/std::log(2.0) <<std::endl;
        //if(print) std::cout<<nb1<<" "<<nb2<<" "<<std::pow(std::pow(2.0,NT(-2)/NT(n)),2)<<std::endl;
        if(print) std::cout<<"\nConstructing the sequence of balls"<<std::endl;

        std::vector<Ball> balls;
        /*
        balls.push_back(Ball(c,std::pow(radius,2)));
        if (print) {
                std::vector<Ball>::iterator bit=balls.end();--bit;
                std::cout<<"ball "<<bit-balls.begin()<<" | "
                         <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
            }
        */
        for(int i=nb1; i<=nb2; ++i){
            balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            if (print) {
                std::vector<Ball>::iterator bit=balls.end();--bit;
                std::cout<<"ball "<<bit-balls.begin()<<" | "<<i
                        <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
            }
        }
        assert(!balls.empty());
        if (print) std::cout<<"---------"<<std::endl;

        // 5. Estimate Vol(P)
        //
        //TODO: std::forward_list<Point> randPoints;
        //std::list<Point> randPoints;
        //randPoints.push_front(p);

        EXACT_NT telescopic_prod=EXACT_NT(1);

        std::vector<Ball>::iterator bit2=balls.end();
        bit2--;

        while(bit2!=balls.begin()){

            //each step starts with some random points in PBLarge stored in list "randPoints"
            //these points have been generated in a previous step

            BallPoly PBLarge(P,*bit2);
            --bit2;
            BallPoly PBSmall(P,*bit2);

            if(print)
                std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()<<") Ball ratio radius="
                        <<PBLarge.second().radius()<<","<<PBSmall.second().radius()<<std::endl;

            // choose a point in PBLarge to be used to generate more rand points
            Point p_gen = *randPoints.begin();

            // num of points in PBSmall and PBLarge
            int nump_PBSmall = 0;
            int nump_PBLarge = randPoints.size();

            if(print) std::cout<<"Points in PBLarge="<<randPoints.size()
                              <<std::endl;

            //keep the points in randPoints that fall in PBSmall
            std::list<Point>::iterator rpit=randPoints.begin();
            while(rpit!=randPoints.end()){
                if (PBSmall.second().is_in(*rpit) == 0){//not in
                    rpit=randPoints.erase(rpit);
                } else {
                    ++nump_PBSmall;
                    ++rpit;
                }
            }

            if(print) std::cout<<"Points in PBSmall="<<randPoints.size()
                              <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                             <<std::endl;

            if(print) std::cout<<"Generate "<<rnum-nump_PBLarge<<  " more "
                              <<std::endl;

            //generate more random points in PBLarge to have "rnum" in total
            rand_point_generator(PBLarge,p_gen,rnum-nump_PBLarge,walk_len,randPoints,PBSmall,nump_PBSmall,var);

            telescopic_prod *= EXACT_NT(rnum)/EXACT_NT(nump_PBSmall);
            if(print) std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                              <<"\ncurrent_vol="<<telescopic_prod
                             <<"\n="<<CGAL::to_double(telescopic_prod)
                            <<"\n--------------------------"<<std::endl;

            //don't continue in pairs of balls that are almost inside P, i.e. ratio ~= 2
            //if(NT(rnum)/NT(nump_PBSmall)>double(1.999)){
            //	break;
            //}
        }
        //if(print) std::cout << "Stopped " << (bit2-balls.begin()) << " balls before Chebychev ball."<< std::endl;
        //telescopic_prod *= std::pow(2,(bit2-balls.begin()));
        if(print) std::cout<<"rand points = "<<rnum<<std::endl;
        if(print) std::cout<<"walk len = "<<walk_len<<std::endl;
        const NT pi = boost::math::constants::pi<NT>();
        //NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0)
        //NT vol = (2*std::pow(pi,n/2.0)*std::pow(radius,n)) / (std::tgamma(n/2.0)*n)

        mpfr_t result,pow,base,exp;
        mpfr_init(result);
        mpfr_init(pow);
        mpfr_init(base);
        mpfr_init(exp);
        mpfr_set_ld(result,2.0,GMP_RNDN);

        mpfr_set_ld(base,pi,GMP_RNDN);
        mpfr_set_ld(exp,n/2.0,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);

        mpfr_set_ld(base,balls[0].radius(),GMP_RNDN);
        mpfr_set_ld(exp,n,GMP_RNDN);
        mpfr_pow(pow, base, exp, GMP_RNDN);
        mpfr_mul(result,result,pow,GMP_RNDN);
        mpfr_div_d(result,result,std::tgamma(n/2.0)*n,GMP_RNDN);
        mpfr_mul_d(result,result,CGAL::to_double(telescopic_prod),GMP_RNDN);

        //std::cout << "mpfr vol=" << mpfr_get_ld(result,GMP_RNDN) << std::endl;
        //
        /*EXACT_NT vol_thread = EXACT_NT(2)
                            * EXACT_NT(std::pow(pi,n/2.0))
                            * EXACT_NT(std::pow(balls[0].radius(),n))
                            / EXACT_NT(EXACT_NT(std::tgamma(n/2.0))*EXACT_NT(n))
                            //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
                            * telescopic_prod;
        */
        //NT vol(0);
        //#pragma omp ordered
        NT vol_thread = mpfr_get_d(result,GMP_RNDN);
        vol += vol_thread;
    }

    // std::cout<<"ROUNDING:"<<round_value<<", "<<CGAL::to_double(round_value*(vol/n_threads)) << ", " <<
    //           CGAL::to_double(round_value*(vol/n_threads)/n*(n+1))<<std::endl;
    //const NT pi = boost::math::constants::pi<NT>();
    //std::cout<<"Cheb:"<<(2*std::pow(pi,n/2.0)*std::pow(radius,n))
    //	                    / (std::tgamma(n/2.0)*n)<<std::endl;
    return round_value*(vol/n_threads);
}
