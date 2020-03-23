// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H




// Pick a random direction as a normilized vector
template <typename RNGType, typename Point, typename NT>
Point get_direction(const unsigned int dim) {

    boost::normal_distribution<> rdist(0,1);
    NT normal = NT(0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    Point p(dim);
    NT* data = p.pointerToData();

    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; ++i) {
        *data = rdist(rng);
        normal += *data * *data;
        data++;
    }

    normal=1.0/std::sqrt(normal);
    p *= normal;

    return p;
}


// Pick a random point from a d-sphere
template <typename RNGType, typename Point, typename NT>
Point get_point_on_Dsphere(const unsigned int dim, const NT &radius){
    Point p = get_direction<RNGType, Point, NT>(dim);
    p = (radius == 0) ? p : radius * p;
    return p;
}


// Pick a random point from a d-ball
template <typename RNGType, typename Point, typename NT>
Point get_point_in_Dsphere(const unsigned int dim, const NT &radius){

    boost::random::uniform_real_distribution<> urdist(0,1);
    NT U;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng2(seed);
    Point p = get_direction<RNGType, Point, NT>(dim);
    U = urdist(rng2);
    U = std::pow(U, 1.0/(NT(dim)));
    p *= (radius*U);
    return p;
}

// ball walk with uniform target distribution
template <typename RNGType, typename Point, typename Polytope, typename NT>
void ball_walk(Point &p,
               Polytope &P,
               const NT &delta)
{
    //typedef typename Parameters::RNGType RNGType;
    Point y = get_point_in_Dsphere<RNGType, Point>(p.dimension(), delta);
    y +=  p;
    if (P.is_in(y)==-1) p = y;
}

// WARNING: USE ONLY WITH BIRKHOFF POLYOPES
// Compute more random points using symmetries of birkhoff polytope
/*
template <class T, class Point, class K, typename  NT>
int birk_sym(T &P, K &randPoints, Point &p) {
    int n=std::sqrt(p.dimension());
    std::vector<int> myints;
    for (unsigned int i=0; i<n; i++){
        myints.push_back(i);
    }

    //The n! possible permutations with n elements
    do {
        std::vector<NT> newpv;

        for (unsigned int j=0; j<p.dimension(); j++){
            int idx = (myints[j/n])*n+1+j%n-1;
            newpv.push_back(p[idx]);
        }

        Point new_p(p.dimension(),newpv.begin(),newpv.end());
        if(P.is_in(new_p) != 0) {
            randPoints.push_back(new_p);
        }
    } while ( std::next_permutation(myints.begin(),myints.end()) );
}*/


// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <typename Polytope, typename PointList, typename Parameters, typename Point>
void boundary_rand_point_generator(Polytope &P,
                                   Point &p,   // a point to start
                                   const unsigned int rnum,
                                   const unsigned int walk_len,
                                   PointList &randPoints,
                                   const Parameters &var)  // constants for volume
{
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    typedef typename Point::Coeff VT;

    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    VT lamdas, Av;

    lamdas.setZero(P.num_of_hyperplanes());
    Av.setZero(P.num_of_hyperplanes());

    unsigned int rand_coord, rand_coord_prev;
    NT kapa, lambda;
    Point p_prev = p, p1(n), p2(n), v(n);
    std::pair <NT, NT> bpair;

    if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else{
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av);
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p += (lambda * v);
    }
    //hit_and_run(p, P, var);

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.cdhr_walk) {

                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
                p_prev = p;
                p1 = p;
                p2 = p;
                p1.set_coord(rand_coord, p[rand_coord] + bpair.first);
                p2.set_coord(rand_coord, p[rand_coord] + bpair.second);
                p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));

            } else {

                v = get_direction<RNGType, Point, NT>(n);
                std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
                p1 = (bpair.first * v) + p;
                p2 = (bpair.second * v) + p;
                lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
                p += (lambda * v);

            }
        }
        randPoints.push_back(p1);
        randPoints.push_back(p2);
    }

}


template <typename Polytope, typename PointList, typename Parameters, typename Point>
void rand_point_generator(Polytope &P,
                         Point &p,   // a point to start
                         const unsigned int rnum,
                         const unsigned int walk_len,
                         PointList &randPoints,
                         const Parameters &var)  // constants for volume
{
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n = var.n;
    RNGType &rng = var.rng; 
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    typedef typename Point::Coeff VT;
    VT lamdas, Av;

    lamdas.setZero(P.num_of_hyperplanes());
    Av.setZero(P.num_of_hyperplanes());

    unsigned int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta, lambda;
    Point p_prev = p, v(n);

    if (var.ball_walk) {
        ball_walk <RNGType> (p, P, ball_rad);
    }else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else if (var.rdhr_walk) {
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av);
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p += (lambda * v);
        //hit_and_run(p, P, var);
    } else {
        billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, true);
    }
    randPoints.push_back(p);

    for (unsigned int i = 1; i <= rnum-1; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk<RNGType>(p, P, ball_rad);
            }else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas);
            } else if (var.rdhr_walk) {
                v = get_direction<RNGType, Point, NT>(n);
                std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
                lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
                p += (lambda * v);
                //hit_and_run(p, P, var);
            } else {
                billiard_walk(P, p, var.diameter, lamdas, Av, lambda,  var);
            }
        }
        randPoints.push_back(p);
    }
 
}



template <typename BallPoly, typename PointList, typename Parameters, typename Point>
void rand_point_generator(BallPoly &PBLarge,
                         Point &p,   // a point to start
                         const unsigned int rnum,
                         const unsigned int walk_len,
                         PointList &randPoints,
                         const BallPoly &PBSmall,
                         unsigned int &nump_PBSmall,
                         const Parameters &var) {  // constants for volume

    typedef typename Point::FT NT;
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    typedef typename Point::Coeff VT;
    VT lamdas, Av;

    lamdas.setZero(PBLarge.num_of_hyperplanes());
    Av.setZero(PBLarge.num_of_hyperplanes());

    unsigned int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta, lambda;
    Point p_prev = p, v(n);

    if (var.ball_walk) {
        ball_walk<RNGType>(p, PBLarge, ball_rad);
    }else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        
        std::pair <NT, NT> bpair = PBLarge.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else if (var.rdhr_walk) {
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = PBLarge.line_intersect(p, v, lamdas, Av);
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p += (lambda * v);
        //hit_and_run(p, PBLarge, var);
    } else {
        billiard_walk(PBLarge, p, var.diameter, lamdas, Av, lambda, var, true);
    }

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk<RNGType>(p, PBLarge, ball_rad);
            }else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, PBLarge, rand_coord, rand_coord_prev, kapa, lamdas);
            } else if (var.rdhr_walk) {
                v = get_direction<RNGType, Point, NT>(n);
                std::pair <NT, NT> bpair = PBLarge.line_intersect(p, v, lamdas, Av, lambda);
                lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
                p += (lambda * v);
                //hit_and_run(p, PBLarge, var);
            } else {
                billiard_walk(PBLarge, p, var.diameter, lamdas, Av, lambda, var);
            }
        }
        if (PBSmall.second().is_in(p) == -1) {//is in
            randPoints.push_back(p);
            ++nump_PBSmall;
        }
    }
}

template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void uniform_first_point(Polytope &P,
                         Point &p,   // a point to start
                         Point &p_prev, // previous point
                         unsigned int &coord_prev, // previous coordinate ray
                         unsigned int walk_len, // number of steps for the random walk
                         VT &lamdas,
                         VT &Av,
                         NT &lambda,
                         const Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n, rand_coord;
    NT kapa, ball_rad = var.delta;
    Point v(n);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    if (var.ball_walk) {
        ball_walk<RNGType>(p, P, ball_rad);
    } else if (var.cdhr_walk) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        coord_prev = rand_coord;
    } else if (var.rdhr_walk) {
        v = get_direction<RNGType, Point, NT>(n);
        std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av);
        lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        p += (lambda * v);
    } else {
        billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, true);
    }
    walk_len--;


    if (var.ball_walk) {
        for (unsigned int j = 0; j < walk_len; j++) {
            ball_walk<RNGType>(p, P, ball_rad);
        }
    } else if (var.cdhr_walk) {
        for (unsigned int j = 0; j < walk_len; j++) {
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, kapa, lamdas);
            coord_prev = rand_coord;
        }
    } else if (var.rdhr_walk) {
        for (unsigned int j = 0; j < walk_len; j++) {
            v = get_direction<RNGType, Point, NT>(n);
            std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
            lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
            p += (lambda * v);
        }
    } else {
        billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var);
    }
}



template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void uniform_next_point(Polytope &P,
                        Point &p,   // a point to start
                        Point &p_prev, // previous point
                        unsigned int &coord_prev, // previous coordinate ray
                        const unsigned int walk_len, // number of steps for the random walk
                        VT &lamdas,
                        VT &Av,
                        NT &lambda,
                        const Parameters &var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n, rand_coord;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    NT ball_rad = var.delta, kapa;
    Point v(n);


    if (var.ball_walk) {
        for (unsigned int j = 0; j < walk_len; j++) ball_walk<RNGType>(p, P, ball_rad);
    } else if (var.rdhr_walk) {
        for (unsigned int j = 0; j < walk_len; j++) {
            v = get_direction<RNGType, Point, NT>(n);
            std::pair <NT, NT> bpair = P.line_intersect(p, v, lamdas, Av, lambda);
            lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
            p += (lambda * v);
        }
    } else if (var.bill_walk) {
        for (unsigned int j = 0; j < walk_len; j++) billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var);
    }else {
        for (unsigned int j = 0; j < walk_len; j++) {
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, kapa, lamdas);
            coord_prev = rand_coord;
        }
    }
}

// ----- HIT AND RUN FUNCTIONS ------------ //
/*
//hit-and-run with random directions and update
template <typename Polytope, typename Point, typename Parameters>
void hit_and_run(Point &p,
                Polytope &P,
                Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    typedef typename Point::FT NT;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point v = get_direction<RNGType, Point, NT>(n);
    std::pair <NT, NT> bpair = P.line_intersect(p, v);
    //NT lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
    p = (( urdist(rng) * (bpair.first - bpair.second) + bpair.second) * v) + p;

}*/


//hit-and-run with orthogonal directions and update
template <typename Polytope, typename Point, typename NT, typename VT>
void hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             Polytope &P,
                             unsigned int rand_coord,
                             unsigned int rand_coord_prev,
                             const NT &kapa,
                             VT &lamdas) {
    std::pair <NT, NT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
}


template <class ConvexBody, class Point, class Parameters, typename NT, typename VT>
void billiard_walk(ConvexBody &P, Point &p, NT diameter, VT &Ar, VT &Av, NT &lambda_prev,
                   Parameters &var, bool first = false) {

    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    NT T = urdist(rng) * diameter;
    const NT dl = 0.995;
    Point v = get_direction<RNGType, Point, NT>(n), p0 = p;
    int it = 0;

    if (first) {

        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av);
        if (T <= pbpair.first) {
            p = (T * v) + p;
            lambda_prev = T;
            return;
        }
        lambda_prev = dl * pbpair.first;
        p += (lambda_prev * v);
        T -= lambda_prev;
        P.compute_reflection(v, p, pbpair.second);
    }

    while (it<10*n) {

        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev);
        if (T <= pbpair.first) {
            p = (T * v) + p;
            lambda_prev = T;
            break;
        }

        lambda_prev = dl * pbpair.first;
        p += (lambda_prev * v);
        T -= lambda_prev;
        P.compute_reflection(v, p, pbpair.second);
        it++;
    }

    if(it == 10*n) p = p0;
}


#endif //RANDOM_SAMPLERS_H
