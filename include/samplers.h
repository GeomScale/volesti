// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H


// Pick a random direction as a normilized vector
Point get_direction(int dim) {
    boost::normal_distribution<> rdist(0,1);
    std::vector<NT> Xs(dim,0);
    NT normal = NT(0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = rdist(rng);
        normal += Xs[i] * Xs[i];
    }
    normal=1.0/std::sqrt(normal);

    for (unsigned int i=0; i<dim; i++) {
        Xs[i] = Xs[i] * normal;
    }
    Point p(dim, Xs.begin(), Xs.end());
    return p;
}


// Pick a random point from a d-sphere
template <typename FT>
Point get_point_on_Dsphere(int dim, FT radius){
    Point p = get_direction(dim);
    p = radius * p;
    return p;
}


// Pick a random point from a d-ball
template <typename FT>
Point get_point_in_Dsphere(int dim, FT radius){
    boost::random::uniform_real_distribution<> urdist(0,1);
    FT U;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng2(seed);
    Point p = get_direction(dim);
    U = urdist(rng2);
    U = std::pow(U, 1.0/(FT(dim)));
    p = (radius*U)*p;
    return p;
}


// WARNING: USE ONLY WITH BIRKHOFF POLYOPES
// Compute more random points using symmetries of birkhoff polytope
template <class T, class K>
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
}


// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <class T, class K>
void rand_point_generator(T &P,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         vars &var)  // constants for volume
{
    int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    std::vector <NT> lamdas(P.num_of_hyperplanes(), NT(0));
    int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta;
    Point p_prev = p;

    if (var.ball_walk) {
        ball_walk(p, P, ball_rad);
    }else if (var.coordinate) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else
        hit_and_run(p, P, var, var);

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk(p, P, ball_rad);
            }else if (var.coordinate) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas, var, var);
            } else
                hit_and_run(p, P, var, var);
        }
        randPoints.push_back(p);
    }
}



template <class T, class K, typename FT>
void rand_point_generator(BallIntersectPolytope<T,FT> &PBLarge,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         BallIntersectPolytope<T,FT> &PBSmall,
                         int &nump_PBSmall,
                         vars &var) {  // constants for volume
    int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    std::vector <FT> lamdas(PBLarge.num_of_hyperplanes(), FT(0));
    int rand_coord, rand_coord_prev;
    NT kapa, ball_rad = var.delta;
    Point p_prev = p;

    if (var.ball_walk) {
        ball_walk(p, PBLarge, ball_rad);
    }else if (var.coordinate) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = PBLarge.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else {
        hit_and_run(p, PBLarge, var, var);
    }

    for (unsigned int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                ball_walk(p, PBLarge, ball_rad);
            }else if (var.coordinate) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, PBLarge, rand_coord, rand_coord_prev, kapa, lamdas, var, var);
            } else
                hit_and_run(p, PBLarge, var, var);
        }
        if (PBSmall.second().is_in(p) == -1) {//is in
            randPoints.push_back(p);
            ++nump_PBSmall;
        }
    }
}

// ----- HIT AND RUN FUNCTIONS ------------ //

//hit-and-run with random directions and update
template <class T>
void hit_and_run(Point &p,
                T &P,
                vars &var,
                vars &var2) {
    int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point l = get_direction(n);
    std::pair <NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
}


//hit-and-run with orthogonal directions and update
template <class T, typename FT>
void hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             T &P,
                             int rand_coord,
                             int rand_coord_prev,
                             FT kapa,
                             std::vector<FT> &lamdas,
                             vars &var,
                             vars &var2) {
    std::pair <FT, FT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
}


// ball walk with uniform target distribution
template <class T, typename FT>
int ball_walk(Point &p,
              T &P,
              FT delta)
{
    Point y = get_point_in_Dsphere(p.dimension(), delta);
    y = y + p;
    if (P.is_in(y)==-1) {
        p = y;
    }
    return 1;
}



#endif //RANDOM_SAMPLERS_H
