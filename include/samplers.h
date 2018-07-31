// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H


template <class P>
class Random_points_on_sphere_d
{
private:
    typedef typename P::FT 	FT;
    int d;
    FT r;
public:
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(std::chrono::system_clock::now().time_since_epoch().count());
    //rn generator;
    //std::default_random_engine generator;   
    //std::normal_distribution<FT> distribution = std::normal_distribution<FT>(FT(0),FT(1));
    boost::normal_distribution<> distribution = boost::normal_distribution<>(0,1.0);
    //typedef std::vector<K> coeff;
    //coeff coeffs;
    
    Random_points_on_sphere_d() {}
    
    Random_points_on_sphere_d(int dim, FT radius) {
        d = dim;
        r = radius;
    }
    
    template <typename GeneratorType>
    P sample_point(GeneratorType generator){
        std::vector<FT> Xs(d,0);
        FT normal=FT(0);
        for (int i=0; i<d; i++){
            Xs[i]=distribution(generator);
            normal+=Xs[i]*Xs[i];
        }
        normal=1.0/std::sqrt(normal);
        
        for (int i=0; i<d; i++){
            Xs[i]=Xs[i]*normal;
        }
        P point(d, Xs.begin(), Xs.end());
        point=point*r;
        
        return point;
    }
    

};

// WARNING: USE ONLY WITH BIRKHOFF POLYOPES
// Compute more random points using symmetries of birkhoff polytope
template <class T, class K>
int birk_sym(T &P,K &randPoints,Point &p){
    int n=std::sqrt(p.dimension());
    std::vector<int> myints;
    for (int i=0; i<n; i++){
        myints.push_back(i);
    }

    //std::cout << "The n! possible permutations with n elements:\n";
    do {
        std::vector<double> newpv;
        for (int i=0; i<n; i++){
            //std::cout << myints[i] << " ";
        }
        //std::cout << std::endl;
        for (int j=0; j<p.dimension(); j++){
            //std::cout << (myints[j/n])*n+1+j%n << " ";
            int idx = (myints[j/n])*n+1+j%n-1;
            //std::cout << idx << " ";
            newpv.push_back(p[idx]);
        }
        //std::cout << "\n";
        Point new_p(p.dimension(),newpv.begin(),newpv.end());
        //std::cout << p << std::endl;
        //std::cout << new_p << "\n" << std::endl;

        //std::cout << P.is_in(new_p) << std::endl;
        if(P.is_in(new_p) != 0){
            //std::cout << "wrong\n";
            randPoints.push_back(new_p);
            //exit(1);
        }
    } while ( std::next_permutation(myints.begin(),myints.end()) );
}


// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <class T, class K>
int rand_point_generator(T &P,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         vars &var)  // constans for volume
{
    int n = var.n;
    bool birk = var.birk;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    std::vector <NT> lamdas(P.num_of_hyperplanes(), NT(0));
    int rand_coord;
    double kapa;
    Point p_prev = p;

    if (var.coordinate) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else
        hit_and_run(p, P, var, var);

    for (int i = 1; i <= rnum; ++i) {

        for (int j = 0; j < walk_len; ++j) {
            int rand_coord_prev = rand_coord;
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            if (var.coordinate) {
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas, var, var);
            } else
                hit_and_run(p, P, var, var);
        }
        randPoints.push_back(p);
    }
}



template <class T, class K>
int rand_point_generator(BallIntersectPolytope<T> &PBLarge,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         BallIntersectPolytope<T> &PBSmall,
                         int &nump_PBSmall,
                         vars &var) {  // constans for volume

    int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    std::vector <NT> lamdas(PBLarge.num_of_hyperplanes(), NT(0));
    int rand_coord;
    double kapa;
    Point p_prev = p;

    if (var.coordinate) {//Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair <NT, NT> bpair = PBLarge.line_intersect_coord(p, rand_coord, lamdas);
        p_prev = p;
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    } else {
        hit_and_run(p, PBLarge, var, var);
    }

    for (int i = 1; i <= rnum; ++i) {
        for (int j = 0; j < walk_len; ++j) {
            int rand_coord_prev = rand_coord;
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            if (var.coordinate)
                hit_and_run_coord_update(p, p_prev, PBLarge, rand_coord, rand_coord_prev, kapa, lamdas, var, var);
            else
                hit_and_run(p, PBLarge, var, var);
        }
        if (PBSmall.second().is_in(p) == -1) {//is in
            randPoints.push_back(p);
            ++nump_PBSmall;
        }
    }
}

// ----- HIT AND RUN FUNCTIONS ------------ //

//hit-and-run with random directions
template <class T>
int hit_and_run(Point &p,
                T &P,
                vars &var,
                vars &var2) {
    int n = var.n;
    //double err = var.err;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);

    Point origin(n);
    Random_points_on_sphere_d<Point> gen(n, 1.0);
    Point l = gen.sample_point(rng);

    std::pair <NT, NT> dbpair = P.line_intersect(p, l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point b1 = (min_plus * l) + p;
    Point b2 = (max_minus * l) + p;
    NT lambda = urdist(rng);
    p = (lambda * b1);
    p = ((1 - lambda) * b2) + p;
    return 1;
}


//hit-and-run with orthogonal directions and update
template <class T>
int hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             T &P,
                             int rand_coord,
                             int rand_coord_prev,
                             double kapa,
                             std::vector<NT> &lamdas,
                             vars &var,
                             vars &var2) {
    std::pair <NT, NT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    p_prev = p;
    p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    return 1;
}



#endif //RANDOM_SAMPLERS_H
