// RandGeom is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// RandGeom is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with RandGeom,
// see <http://www.gnu.org/licenses/>.
//
// Developer: Vissarion Fisikopoulos


#include <CGAL/point_generators_d.h>
//#include <CGAL/Filtered_kernel_d.h>
//#include <CGAL/Triangulation.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/algorithm.h>
//#include <CGAL/Random.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <forward_list>
#include <list>
#include <bitset>
#include <random>
#include <chrono>       // std::chrono::system_clock
#include <functional>
#include <algorithm>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "boost/math/constants/constants.hpp"
#include "boost/dynamic_bitset.hpp"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>
#include <vector>
#include <iostream>

#include <Eigen/Eigen>
//#include <Eigen/Cholesky>

//#include <CGAL/Extreme_points_d.h>
//#include <CGAL/Extreme_points_traits_d.h>

//#include <gmpxx.h>
//typedef mpq_class NT;
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpq                  EXACT_NT;
typedef double                      NT;
//typedef CGAL::Gmpz                NT;

typedef CGAL::Cartesian_d<NT> 	      Kernel;
//typedef CGAL::Triangulation<Kernel> T;
typedef Kernel::Point_d								Point;
typedef Kernel::Vector_d							Vector;
typedef Kernel::Line_d								Line;
typedef Kernel::Hyperplane_d					Hyperplane;
typedef Kernel::Ray_d								Ray;
typedef Kernel::Direction_d						Direction;
//typedef Kernel::Sphere_d						Ball;
typedef CGAL::Approximate_min_ellipsoid_d_traits_d<Kernel, EXACT_NT> Traits;
//typedef Traits::Point                                          Point;
//typedef std::vector<Point>                                     Point_list;
typedef CGAL::Approximate_min_ellipsoid_d<Traits>              AME;

// define random generator
//typedef boost::mt11213b RNGType; ///< mersenne twister generator
typedef boost::mt19937 RNGType; ///< mersenne twister generator
//typedef boost::lagged_fibonacci607 RNGType;
//typedef boost::hellekalek1995 RNGType;
//typedef boost::rand48 RNGType;
//typedef boost::minstd_rand RNGType;

typedef boost::variate_generator< RNGType, boost::normal_distribution<> >  generator;
//typedef boost::variate_generator< RNGType, boost::exponential_distribution<> >  generator;

//structs with variables and random generators
struct vars {
public:
    vars( int m,
          int n,
          int walk_steps,
          int n_threads,
          const double err,
          const double err_opt,
          const int lw,
          double up,
          const int L,
          RNGType &rng,
          generator
          &get_snd_rand,
          boost::random::uniform_real_distribution<> urdist,
          boost::random::uniform_real_distribution<> urdist1,
          bool verbose,
          bool rand_only,
          bool round,
          bool NN,
          bool birk,
          bool coordinate,
          bool use_jl=true,
          bool epsilon=0.1
        ) :
        m(m), n(n), walk_steps(walk_steps), n_threads(n_threads), err(err), err_opt(err_opt),
        lw(lw), up(up), L(L), rng(rng), get_snd_rand(get_snd_rand),
        urdist(urdist), urdist1(urdist1) , verbose(verbose), rand_only(rand_only), round(round),
        NN(NN),birk(birk),coordinate(coordinate),epsilon(epsilon),use_jl(use_jl) {};

    int m;
    int n;
    double epsilon;
    bool use_jl;
    int walk_steps;
    int n_threads;
    const double err;
    const double err_opt;
    const int lw;
    double up;
    const int L;
    RNGType &rng;
    generator
    &get_snd_rand;
    boost::random::uniform_real_distribution<> urdist;
    boost::random::uniform_real_distribution<> urdist1;
    bool verbose;
    bool rand_only;
    bool round;
    bool NN;
    bool birk;
    bool coordinate;
};

// define extreme points
//typedef CGAL::Extreme_points_traits_d<Point>   EP_Traits_d;

template <class T>
int optimization(T &KK,vars var,Point &fp,Vector &w);
template <class T>
int opt_interior(T &K,vars &var,Point &opt,Vector &w);

#include <polytopes.h>
#include <ballintersectpolytope.h>
//#include <opt_rand.h>
//#include <oracles.h>
//#include <../../vol2/include/random_samplers.h>
#include <random_samplers.h>
#include <rounding.h>
#include <misc.h>
#include <linear_extensions.h>

/////////////////////////////////////////////////////////
// VOLUME
// randomized approximate volume computation
/*************************************************
/* VOLUME with random DIRECTIONS hit and run     */
// We assume that the polytope P is properly sandwitched
// The sandwitching:
// r is the radius of the smallest ball
// d is the radius of the largest
template <class T>
NT volume0(T &P,
           vars &var,  // constans for volume
           vars &var2, // constants for optimization in case of MinkSums
           NT r,
           NT d) {
    typedef BallIntersectPolytope<T>        BallPoly;

    bool print = true;
    int n = var.n;
    int rnum = var.m;
    int walk_len = var.walk_steps;
    const double err = var.err;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

    //The number of balls we construct
    int nb = std::ceil(n * (std::log(d)/std::log(2.0)));
    //std::cout<<"nb="<<nb<<", d="<<d<<std::endl;
    //std::pow(std::log(n),2)

    // Construct the sequence of balls
    std::vector<NT> coords(n,0);
    Point p0(n,coords.begin(),coords.end());
    std::vector<Ball> balls;
    for(int i=0; i<=nb; ++i) {
        balls.push_back(Ball(p0,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
        if (print) std::cout<<"ball"<<i<<"="<<balls[i].center()<<" "<<balls[i].radius()<<std::endl;
    }
    assert(!balls.empty());
    if (print) std::cout<<"---------"<<std::endl;


    std::vector<int> telescopic_prod(nb,0);
    for(int i=1; i<=rnum; ++i) { //generate rnum rand points
        //start with a u.d.r point in the smallest ball
        //radius=1, center=Origin()
        std::vector<NT> coords(n,0);
        Point p(n,coords.begin(),coords.end());
        BallPoly PBold(P,balls[0]);
        //
        //hit_and_run(p,PBold.second(),var,var2);
        CGAL::Random_points_in_ball_d<Point> gen (n, 1.0);
        p = *gen;
        //std::cout<<p<<std::endl;
        //std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
        //std::cout<<balls[0].is_in(p)<<std::endl;

        //exit(0);
        std::vector<Ball>::iterator bit=balls.begin();
        std::vector<int>::iterator prod_it=telescopic_prod.begin();
        ++bit;
        for(; bit!=balls.end(); ++bit, ++prod_it) {
            // generate a random point in bit intersection with P
            BallPoly PB(P,*bit);

            for(int j=0; j<walk_len; ++j) {
                hit_and_run(p,PB,var,var2);
                //std::cout<<"h-n-r:"<<p<<std::endl;
            }
            //Not need to test for PBold membership. Just check if inside Ball
            //if (Sep_Oracle(PBold,p,var2).get_is_in()){
            if (PBold.second().is_in(p)) {
                //std::cout<<p<<" IN ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
                ++(*prod_it);
            } else {
                ;
                //std::cout<<p<<":"<<(p-CGAL::Origin()).squared_length()
                //<<" OUT ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
            }
            PBold=PB;
        }
        if (print) std::cout<<"\n\ngenerated random point..."<<i<<"/"<<rnum<<" ";
        const NT pi = boost::math::constants::pi<NT>();
        NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0);
        for(std::vector<int>::iterator prod_it=telescopic_prod.begin();
                prod_it!=telescopic_prod.end(); ++prod_it) {
            vol *= NT(i)/NT(*prod_it);
        }
        if (print) std::cout<<"current vol estimation= "<<vol<<std::endl;
        if (print) std::cout<<"walklen="<<walk_len<<std::endl;
        if (print) std::cout<<"rnum="<<rnum<<std::endl;
        //for(prod_it=telescopic_prod.begin(); prod_it!=telescopic_prod.end(); ++prod_it)
        //	std::cout<<(*prod_it)<<" ";
        //std::cout<<std::endl;
    }
    const NT pi = boost::math::constants::pi<NT>();
    NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0);
    //NT vol=1;
    if (print) std::cout<<"vol(K_0)="<<vol<<" ";
    for(std::vector<int>::iterator prod_it=telescopic_prod.begin();
            prod_it!=telescopic_prod.end(); ++prod_it) {
        vol *= NT(rnum)/NT(*prod_it);
        if (print) std::cout<<NT(rnum)<<"/" << NT(*prod_it)<<"="<<NT(rnum)/NT(*prod_it)<<"\n";
    }
    return vol;
}


/*************************************************
/* VOLUME with random COORDINATES hit and run    */
// We assume that the polytope P is properly sandwitched
// The sandwitching:
// r is the radius of the smallest ball
// d is the radius of the largest
template <class T>
NT volume1(T &P,
           vars &var,  // constans for volume
           vars &var2, // constants for optimization in case of MinkSums
           NT r,
           NT d) {
    typedef BallIntersectPolytope<T>        BallPoly;

    bool print = var.verbose;
    int n = var.n;
    int rnum = var.m;
    int walk_len = var.walk_steps;
    const double err = var.err;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    boost::random::uniform_int_distribution<> uidist(0,n-1);
    //boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

    //The number of balls we construct
    int nb = std::ceil(n * (std::log(d)/std::log(2.0)));
    //std::cout<<"nb="<<nb<<", d="<<d<<std::endl;
    //std::pow(std::log(n),2)

    // Construct the sequence of balls
    std::vector<NT> coords(n,0);
    Point p0(n,coords.begin(),coords.end());
    std::vector<Ball> balls;
    for(int i=0; i<=nb; ++i) {
        balls.push_back(Ball(p0,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
        if (print) std::cout<<"ball"<<i<<"="<<balls[i].center()<<" "<<balls[i].radius()<<std::endl;
    }
    assert(!balls.empty());
    if (print) std::cout<<"---------"<<std::endl;

    std::vector<int> telescopic_prod(nb,0);

    //#pragma omp parallel for ordered schedule(dynamic)
    for(int i=1; i<=rnum; ++i) { //generate rnum rand points
        //start with a u.d.r point in the smallest ball
        //radius=1, center=Origin()
        std::vector<NT> coords(n,0);
        Point p(n,coords.begin(),coords.end());
        BallPoly PBold(P,balls[0]);
        //
        //hit_and_run(p,PBold.second(),var,var2);
        CGAL::Random_points_in_ball_d<Point> gen (n, NT(1.0));
        p = *gen;
        //std::cout<<p<<std::endl;
        //std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
        //std::cout<<balls[0].is_in(p)<<std::endl;

        std::vector<Ball>::iterator bit=balls.begin();
        std::vector<int>::iterator prod_it=telescopic_prod.begin();
        ++bit;
        for(; bit!=balls.end(); ++bit, ++prod_it) {
            // generate a random point in bit intersection with P
            BallPoly PB(P,*bit);

            std::vector<NT> lamdas(P.size(),NT(0));
            int rand_coord = uidist(rng);
            double kapa = urdist(rng);
            Point p_prev=p;
            hit_and_run_coord_update(p,p_prev,PB,rand_coord,rand_coord,kapa,lamdas,var,var2,true);

            for(int j=0; j<walk_len; ++j) {
                int rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p,p_prev,PB,rand_coord,rand_coord_prev,kapa,lamdas,var,var2,false);
            }

            //Not need to test for PBold membership. Just check if inside Ball
            //if (Sep_Oracle(PBold,p,var2).get_is_in()){
            if (PBold.second().is_in(p)) {
                //std::cout<<p<<" IN ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
                ++(*prod_it);
            } else {
                ;
                //std::cout<<p<<":"<<(p-CGAL::Origin()).squared_length()
                //<<" OUT ball: "<<PBold.second().center()<<PBold.second().radius()<<std::endl;
            }
            PBold=PB;
        }
        if (print) std::cout<<"\n\ngenerated random point..."<<i<<"/"<<rnum<<" ";
        const NT pi = boost::math::constants::pi<NT>();
        NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0);
        for(std::vector<int>::iterator prod_it=telescopic_prod.begin();
                prod_it!=telescopic_prod.end(); ++prod_it) {
            vol *= NT(i)/NT(*prod_it);
        }
        if (print) std::cout<<"current vol estimation= "<<vol<<std::endl;
        if (print) std::cout<<"walklen="<<walk_len<<std::endl;
        //if (print) std::cout<<"rnum="<<rnum<<std::endl;
        //for(prod_it=telescopic_prod.begin(); prod_it!=telescopic_prod.end(); ++prod_it)
        //	std::cout<<(*prod_it)<<" ";
        //std::cout<<std::endl;
    }
    const NT pi = boost::math::constants::pi<NT>();
    NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0);
    //NT vol=1;
    if (print) std::cout<<"vol(K_0)="<<vol<<" ";
    for(std::vector<int>::iterator prod_it=telescopic_prod.begin();
            prod_it!=telescopic_prod.end(); ++prod_it) {
        vol *= NT(rnum)/NT(*prod_it);
        if (print) std::cout<<NT(rnum)<<"/" << NT(*prod_it)<<"="<<NT(rnum)/NT(*prod_it)<<"\n";
    }
    return vol;
}

/*************************************************
/* VOLUME with random COORDINATES hit and run
 * Here we reuse the random points we generate   */
// We assume that the polytope P is properly sandwitched
// The sandwitching:
// r is the radius of the smallest ball
// d is the radius of the largest
template <class T>
NT volume1_reuse(T &P,
                 vars &var,  // constans for volume
                 vars &var2, // constants for optimization in case of MinkSums
                 NT r,
                 NT d) {
    typedef BallIntersectPolytope<T>        BallPoly;

    bool print = var.verbose;
    int n = var.n;
    int rnum = var.m;
    int walk_len = var.walk_steps;
    const double err = var.err;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    boost::random::uniform_int_distribution<> uidist(0,n-1);
    //boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

    //The number of balls we construct
    int nb = std::ceil(n * (std::log(d)/std::log(2.0)));

    // Construct the sequence of balls
    std::vector<NT> coords0(n,0);
    Point p0(n,coords0.begin(),coords0.end());
    std::vector<Ball> balls;
    for(int i=0; i<=nb; ++i) {
        balls.push_back(Ball(p0,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
        if (print) std::cout<<"ball"<<i<<"="<<balls[i].center()<<" "<<balls[i].radius()<<std::endl;
    }
    assert(!balls.empty());
    if (print) std::cout<<"---------"<<std::endl;

    //Generate the first random point in P
    //1. start with a u.d.r point in the smallest ball
    std::vector<NT> coords(n,0);
    Point p(n,coords.begin(),coords.end());
    BallPoly PBold(P,balls[0]);
    CGAL::Random_points_in_ball_d<Point> gen (n, NT(1.0));
    p = *gen;
    //2. use p to generate the next rand point in the next ball
    //   until we reach P
    std::vector<Ball>::iterator bit=balls.begin();
    //std::vector<int>::iterator prod_it=telescopic_prod.begin();
    ++bit;
    for(; bit!=balls.end(); ++bit) {
        // generate a random point in bit intersection with P
        BallPoly PB(P,*bit);

        std::vector<NT> lamdas(P.size(),NT(0));
        int rand_coord = uidist(rng);
        double kapa = urdist(rng);
        Point p_prev=p;
        hit_and_run_coord_update(p,p_prev,PB,rand_coord,rand_coord,kapa,lamdas,var,var2,true);

        for(int j=0; j<walk_len; ++j) {
            int rand_coord_prev = rand_coord;
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            hit_and_run_coord_update(p,p_prev,PB,rand_coord,rand_coord_prev,kapa,lamdas,var,var2,false);
        }
        //Not need to test for PBold membership. Just check if inside Ball
        //if (PBold.second().is_in(p))
        //  ++(*prod_it);
        PBold=PB;
    }

    //Now that p is a random point in P
    //use it to generate random points in P
    //TODO: std::forward_list<Point> randPoints;
    std::list<Point> randPoints;
    randPoints.push_front(p);

    NT telescopic_prod=NT(1);

    std::vector<Ball>::iterator bit2=balls.end();
    bit2--;


    while(bit2!=balls.begin()) {

        //each step starts with some random points in PBLarge stored in list "randPoints"
        //these points have been generated in a previous step

        BallPoly PBLarge(P,*bit2);
        --bit2;
        BallPoly PBSmall(P,*bit2);

        if (print) std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()<<") Ball ratio radius="
                                <<PBLarge.second().radius()<<","<<PBSmall.second().radius()<<std::endl;

        // choose a point in PBLarge to be used to generate more rand points
        Point p_gen = *randPoints.begin();

        // num of points in PBSmall and PBLarge
        int nump_PBSmall = 0;
        int nump_PBLarge = randPoints.size();

        if (print) std::cout<<"Points in PBLarge="<<randPoints.size()
                                <<std::endl;

        //keep the points in randPoints that fall in PBSmall
        std::list<Point>::iterator rpit=randPoints.begin();
        while(rpit!=randPoints.end()) {
            if (PBSmall.second().is_in(*rpit) == 0) { //not in
                rpit=randPoints.erase(rpit);
            } else {
                ++nump_PBSmall;
                ++rpit;
            }
        }

        if (print) std::cout<<"Points in PBSmall="<<randPoints.size()
                                <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                                <<std::endl;

        if (print) std::cout<<"Generate "<<rnum-nump_PBLarge<<  " more "
                                <<std::endl;

        //generate more random points in PBLarge to have "rnum" in total
        for(int i=1; i<=rnum - nump_PBLarge; ++i) {
            std::vector<NT> lamdas(P.size(),NT(0));
            int rand_coord = uidist(rng);
            double kapa = urdist(rng);
            Point p_gen_prev = p_gen;
            hit_and_run_coord_update(p_gen,p_gen_prev,PBLarge,rand_coord,rand_coord,kapa,lamdas,var,var2,true);
            for(int j=0; j<walk_len; ++j) {
                int rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p_gen,p_gen_prev,PBLarge,rand_coord,rand_coord_prev,kapa,lamdas,var,var2,false);
            }
            // count and store in randPoints the points fall in PBSmall
            if (PBSmall.second().is_in(p_gen) == -1) { //is in
                randPoints.push_back(p_gen);
                ++nump_PBSmall;
            }
        }
        telescopic_prod *= NT(rnum)/NT(nump_PBSmall);
        if (print) std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                                <<"\n--------------------------"<<std::endl;
    }
    if (print) std::cout<<"rand points = "<<rnum<<std::endl;
    if (print) std::cout<<"walk len = "<<walk_len<<std::endl;
    const NT pi = boost::math::constants::pi<NT>();
    NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0)
             //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
             * telescopic_prod;
    //NT vol(0);
    return vol;
}

/*************************************************
/* VOLUME with random COORDINATES hit and run
 * THIS IS A TEST FUNCTION
 *
 * Here we reuse the random points we generate   */
// We assume that the polytope P is properly sandwitched
// The sandwitching:
// r is the radius of the smallest ball
// d is the radius of the largest
template <class T>
NT volume1_reuse_test(T &P,
                      vars &var,  // constans for volume
                      vars &var2 // constants for optimization in case of MinkSums
                     ) {
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
    if(round) {
        round_value = rounding(P,var,var2);
    }

    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    //1. Compute the Chebychev ball (largest inscribed ball) with center and radius
    double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    if(print) std::cout<<"\nComputing the Chebychev center..."<<std::endl;
    Point c;       //center
    double radius;
    P.chebyshev_center(c,radius);
    //radius=std::sqrt(radius);
    if(print) std::cout<<"Chebychev center= "<<c<<"\nradius="<<radius<<std::endl;
    double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
    //Chebtime = tstop - tstart;
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(print) std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;

    rnum=rnum/n_threads;
    //NT vol=0;

    // Perform the procedure for a number of threads and then take the average
    //#pragma omp for ordered schedule(dynamic)
    //for(int t=0; t<n_threads; t++){
    // 2. Generate the first random point in P
    // Perform random walk on random point in the Chebychev ball
    if(print) std::cout<<"\nGenerate the first random point in P"<<std::endl;
    CGAL::Random_points_in_ball_d<Point> gen (n, radius);
    Point p = *gen;
    p = p + (c-CGAL::Origin());
    std::list<Point> randPoints; //ds for storing rand points
    //use a large walk length e.g. 1000
    //rand_point_generator(P, p, 1, 1000, randPoints, var);

    std::vector<NT> lamdas(P.num_of_hyperplanes(),NT(0));
    int rand_coord = uidist(rng);
    double kapa = urdist(rng);
    Point p_prev = p;
    hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord,kapa,lamdas,var,var,true);
    for(int j=0; j<1000; ++j) {
        int rand_coord_prev = rand_coord;
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        //hit_and_run(p,P,var,var2);
        hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord_prev,kapa,lamdas,var,var,false);
    }
    randPoints.push_back(p);


    if (print) std::cout<<"First random point: "<<p<<std::endl;

    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
    // 3. Sample "rnum" points from P
    if(print) std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
    //randPoints.push_front(p);
    //rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);

    for(int i=1; i<=rnum-1; ++i) {
        std::vector<NT> lamdas(P.num_of_hyperplanes(),NT(0));
        int rand_coord = uidist(rng);
        double kapa = urdist(rng);
        Point p_prev = p;
        hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord,kapa,lamdas,var,var,true);

        for(int j=0; j<walk_len; ++j) {
            int rand_coord_prev = rand_coord;
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            //hit_and_run(p,P,var,var2);
            hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord_prev,kapa,lamdas,var,var,false);
        }
        randPoints.push_back(p);
    }

    double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(print) std::cout << "First random points construction time = " << tstop2 - tstart2 << std::endl;
    //if(rand_only) return -1;

    // 4.  Construct the sequence of balls
    // 4a. compute the radius of the largest ball
    double current_dist, max_dist=NT(0);
    for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit) {
        current_dist=(*pit-c).squared_length();
        if(current_dist>max_dist) {
            max_dist=current_dist;
        }
    }
    max_dist=std::sqrt(max_dist);
    std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist<<std::endl;

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
    for(int i=nb1; i<=nb2; ++i) {
        balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
        if (print) {
            std::vector<Ball>::iterator bit=balls.end();
            --bit;
            std::cout<<"ball "<<bit-balls.begin()<<" | "<<i
                     <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
        }
    }
    assert(!balls.empty());
    if (print) std::cout<<"---------"<<std::endl;

    //Now that p is a random point in P
    //use it to generate random points in P
    //TODO: std::forward_list<Point> randPoints;
    //std::list<Point> randPoints;
    //randPoints.push_front(p);

    NT telescopic_prod=NT(1);

    std::vector<Ball>::iterator bit2=balls.end();
    bit2--;


    while(bit2!=balls.begin()) {

        //each step starts with some random points in PBLarge stored in list "randPoints"
        //these points have been generated in a previous step

        BallPoly PBLarge(P,*bit2);
        --bit2;
        BallPoly PBSmall(P,*bit2);

        if (print) std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()<<") Ball ratio radius="
                                <<PBLarge.second().radius()<<","<<PBSmall.second().radius()<<std::endl;

        // choose a point in PBLarge to be used to generate more rand points
        Point p_gen = *randPoints.begin();

        // num of points in PBSmall and PBLarge
        int nump_PBSmall = 0;
        int nump_PBLarge = randPoints.size();

        if (print) std::cout<<"Points in PBLarge="<<randPoints.size()
                                <<std::endl;

        //keep the points in randPoints that fall in PBSmall
        std::list<Point>::iterator rpit=randPoints.begin();
        while(rpit!=randPoints.end()) {
            if (PBSmall.second().is_in(*rpit) == 0) { //not in
                rpit=randPoints.erase(rpit);
            } else {
                ++nump_PBSmall;
                ++rpit;
            }
        }

        if (print) std::cout<<"Points in PBSmall="<<randPoints.size()
                                <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                                <<std::endl;

        if (print) std::cout<<"Generate "<<rnum-nump_PBLarge<<  " more "
                                <<std::endl;

        //generate more random points in PBLarge to have "rnum" in total
        for(int i=1; i<=rnum - nump_PBLarge; ++i) {
            std::vector<NT> lamdas(P.num_of_hyperplanes(),NT(0));
            int rand_coord = uidist(rng);
            double kapa = urdist(rng);
            Point p_gen_prev = p_gen;
            hit_and_run_coord_update(p_gen,p_gen_prev,PBLarge,rand_coord,rand_coord,kapa,lamdas,var,var2,true);
            for(int j=0; j<walk_len; ++j) {
                int rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p_gen,p_gen_prev,PBLarge,rand_coord,rand_coord_prev,kapa,lamdas,var,var2,false);
            }
            // count and store in randPoints the points fall in PBSmall
            if (PBSmall.second().is_in(p_gen) == -1) { //is in
                randPoints.push_back(p_gen);
                ++nump_PBSmall;
            }
        }
        telescopic_prod *= NT(rnum)/NT(nump_PBSmall);
        if (print) std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                                <<"\n--------------------------"<<std::endl;
    }
    if (print) std::cout<<"rand points = "<<rnum<<std::endl;
    if (print) std::cout<<"walk len = "<<walk_len<<std::endl;
    const NT pi = boost::math::constants::pi<NT>();

    NT vol = (2*std::pow(pi,n/2.0)*std::pow(balls[0].radius(),n))
             / (std::tgamma(n/2.0)*n)
             //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
             * telescopic_prod;

    //NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0)
    //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
    //       * telescopic_prod;
    //NT vol(0);
    return vol;
}






/*************************************************
/* VOLUME with random COORDINATES hit and run
 * Here we reuse the random points we generate
 * We also use Chebychev ball and sampling for
 * sandwitching
 *   */
template <class T>
NT volume1_reuse2(T &P,
                  vars &var,  // constans for volume
                  vars &var2, // constants for optimization in case of MinkSums
                  double &Chebtime) {
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
    if(round) {
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
    for(int t=0; t<n_threads; t++) {
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
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit) {
            current_dist=(*pit-c).squared_length();
            if(current_dist>max_dist) {
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
        for(int i=nb1; i<=nb2; ++i) {
            balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            if (print) {
                std::vector<Ball>::iterator bit=balls.end();
                --bit;
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

        while(bit2!=balls.begin()) {

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
            while(rpit!=randPoints.end()) {
                if (PBSmall.second().is_in(*rpit) == 0) { //not in
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

        mpfr_clear(result);
        mpfr_clear(pow);
        mpfr_clear(base);
        mpfr_clear(exp);
    }

    // std::cout<<"ROUNDING:"<<round_value<<", "<<CGAL::to_double(round_value*(vol/n_threads)) << ", " <<
    //           CGAL::to_double(round_value*(vol/n_threads)/n*(n+1))<<std::endl;
    //const NT pi = boost::math::constants::pi<NT>();
    //std::cout<<"Cheb:"<<(2*std::pow(pi,n/2.0)*std::pow(radius,n))
    //	                    / (std::tgamma(n/2.0)*n)<<std::endl;
    return round_value*(vol/n_threads);
}



/*************************************************
/* VOLUME with random COORDINATES hit and run
 * Here we reuse the random points we generate
 * We also use Chebychev ball and sampling for
 * sandwitching
 * We also we estimators for the length of the
 * random walk.
 *   */
/*
template <class T>
EXACT_NT volume1_reuse_estimete_walk(T &P,
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
  //radius=std::sqrt(radius);
  if(print) std::cout<<"Chebychev center= "<<c<<"\nradius="<<radius<<std::endl;
  double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
  Chebtime = tstop - tstart;
  double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
	if(print) std::cout << "Chebychev time = " << tstop1 - tstart1 << std::endl;

	rnum=rnum/n_threads;
	EXACT_NT vol=0;

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
		rand_point_generator(P, p, 1, 1000, randPoints, var);
		if (print) std::cout<<"First random point: "<<p<<std::endl;

		double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
		// 3. Sample "rnum" points from P
		if(print) std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
		//randPoints.push_front(p);
		//rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);
		rand_point_generator_with_walk_estimator(P, rnum-1, randPoints, var);
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
		std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist<<std::endl;

		//
		// 4b. Number of balls
		int nb1 = n * (std::log(radius)/std::log(2.0));
		//int nb1 = std::sqrt(n);// * (std::log(radius)/std::log(2.0));
		int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));
		//int nb1 = n * (std::log(radius)/std::log(2.0));
		//int nb2 = n * (std::log(max_dist)/std::log(2.0));
		//std::cout<<n* std::log(radius)/std::log(2.0) <<std::endl;
		//std::cout<<n* std::log(max_dist)/std::log(2.0) <<std::endl;
		//if(print) std::cout<<nb1<<" "<<nb2<<" "<<std::pow(std::pow(2.0,NT(-2)/NT(n)),2)<<std::endl;
		if(print) std::cout<<"\nConstructing the sequence of balls"<<std::endl;

		std::vector<Ball> balls;

		//balls.push_back(Ball(c,std::pow(radius,2)));
		//if (print) {
		//		std::vector<Ball>::iterator bit=balls.end();--bit;
		//		std::cout<<"ball "<<bit-balls.begin()<<" | "
		//		         <<" center="<<bit->center()<<" radius="<<bit->radius()<<std::endl;
		//	}

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
		EXACT_NT vol_thread = EXACT_NT(2*std::pow(pi,n/2.0)*std::pow(balls[0].radius(),n))
		                    / EXACT_NT(std::tgamma(n/2.0)*n)
		                    //* (std::pow(NT(rnum),balls.size()-1) / telescopic_prod_nom );
		                    * telescopic_prod;
		//NT vol(0);
		//#pragma omp ordered
		vol += vol_thread;
	}

	// std::cout<<"ROUNDING:"<<round_value<<", "<<CGAL::to_double(round_value*(vol/n_threads)) << ", " <<
	//           CGAL::to_double(round_value*(vol/n_threads)/n*(n+1))<<std::endl;
	const NT pi = boost::math::constants::pi<NT>();
	//std::cout<<"Cheb:"<<(2*std::pow(pi,n/2.0)*std::pow(radius,n))
	//	                    / (std::tgamma(n/2.0)*n)<<std::endl;
	return round_value*(vol/n_threads);
}
*/


// VOLUME with multipoint random walk
template <class T>
NT volume2(T &P,
           vars &var) {
    typedef BallIntersectPolytope<T>        BallPoly;

    int n = var.n;
    int rnum = var.m;
    int walk_len = var.walk_steps;
    const double err = var.err;
    RNGType &rng = var.rng;
    generator
    get_snd_rand = var.get_snd_rand;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    boost::random::uniform_real_distribution<> urdist1 = var.urdist1;

    // The sandwitching
    // r is the radius of the smallest ball
    // d is the radius of the largest
    std::vector<NT> coords_apex(n,1);
    Vector p_apex(n,coords_apex.begin(),coords_apex.end());
    const NT r=1, d=std::sqrt(p_apex.squared_length());
    const int nb = std::ceil(n * (std::log(d)/std::log(2.0)));
    std::cout<<"nb="<<nb<<", d="<<d<<std::endl;
    //std::pow(std::log(n),2)

    // Construct the sequence of balls
    std::vector<NT> coords(n,0);
    Point p0(n,coords.begin(),coords.end());
    std::vector<Ball> balls;
    for(int i=0; i<=nb; ++i) {
        balls.push_back(Ball(p0,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
        std::cout<<"ball"<<i<<"="<<balls[i].center()<<" "<<balls[i].radius()<<std::endl;
    }
    assert(!balls.empty());
    std::cout<<"---------"<<std::endl;



    std::vector<int> telescopic_prod(nb,0);
    //vector to store the random points
    std::vector<Point> V;
    BallPoly PBold(P,balls[0]);
    for(int i=0; i<rnum; ++i) {
        // generate rnum rand points
        // in the smallest ball i.e. radius=1, center=Origin()
        std::vector<NT> coords(n,0);
        Point p(n,coords.begin(),coords.end());
        hit_and_run(p,PBold,var);
        V.push_back(p);
        //std::cout<<p<<std::endl;
    }
    //std::cout<<p<<std::endl;
    //std::cout<<Sep_Oracle(PBold,p).get_is_in()<<std::endl;
    //std::cout<<balls[0].is_in(p)<<std::endl;
    //exit(0);
    std::vector<Ball>::iterator bit=balls.begin();
    std::vector<int>::iterator prod_it=telescopic_prod.begin();
    ++bit;
    for(; bit!=balls.end(); ++bit, ++prod_it) {
        // generate a random point in bit (intersection) P
        BallPoly PB(P,*bit);
        std::cout<<"walking..."<<walk_len<<"steps"<<std::endl;
        var.m = V.size();
        multipoint_random_walk(PB,V,var);

        for(int j=0; j<V.size(); ++j) {

            if (Sep_Oracle(PBold,V[j],var).get_is_in()) {
                //std::cout<<p<<" IN ball: "<<PBold.second.center()<<PBold.second.radius()<<std::endl;
                ++(*prod_it);
            }//else{
            //std::cout<<p<<":"<<(p-CGAL::Origin()).squared_length()
            //<<" OUT ball: "<<PBold.second.center()<<PBold.second.radius()<<std::endl;
        }
        PBold=PB;
        //for(prod_it=telescopic_prod.begin(); prod_it!=telescopic_prod.end(); ++prod_it)
        //	std::cout<<(*prod_it)<<" ";
        //std::cout<<std::endl;
    }
    const NT pi = boost::math::constants::pi<NT>();
    NT vol = std::pow(pi,n/2.0)/std::tgamma(1+n/2.0);
    //NT vol=1;
    std::cout<<"vol(K_0)="<<vol<<" ";
    for(std::vector<int>::iterator prod_it=telescopic_prod.begin();
            prod_it!=telescopic_prod.end(); ++prod_it) {
        vol *= NT(rnum)/NT(*prod_it);
        std::cout<<NT(rnum)<<"/" << NT(*prod_it)<<"="<<NT(rnum)/NT(*prod_it)<<"\n";
    }
    std::cout<<std::endl;
    std::cout<<"volume = "<<vol<<std::endl;
    std::cout<<"exact volume = "<<std::pow(2,n)<<std::endl;
    return vol;
}

