// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef GAUSSIAN_SAMPLERS_H
#define GAUSSIAN_SAMPLERS_H

int rand_exp_range(NT l, NT u, NT a_i, NT &dis, vars &var){
    NT r, r_val, fn;
    if(a_i>std::pow(10,-8.0) && std::abs(u-l)>=2.0/std::sqrt(2.0*a_i)){
        boost::normal_distribution<> rdist(-1,1);
        while(true){
            r = rdist(var.rng);
            r = r/std::sqrt(2.0*a_i);
            if(r>=l && r<=u){
                break;
            }
        }
        dis=r;

    }else{
        boost::random::uniform_real_distribution<> urdist(0,1);
        NT M=1.0;
        while(true){
            r=urdist(var.rng);
            dis = (1.0-r)*l+r*u;
            r_val = M*urdist(var.rng);
            fn = std::exp(-a_i*dis*dis);
            if(r_val<fn){
                break;
            }
        }
    }
    return 1;
}


int rand_exp_range_coord(NT l, NT u, NT a_i, NT &dis, vars &var){
    NT r, r_val, fn;
    if(a_i>std::pow(10,-8.0) && std::abs(u-l)>=2.0/std::sqrt(2.0*a_i)){
        boost::normal_distribution<> rdist(-1,1);
        while(true){
            r = rdist(var.rng);
            r = r/std::sqrt(2.0*a_i);
            if(r>=l && r<=u){
                break;
            }
        }
        dis=r;

    }else{
        boost::random::uniform_real_distribution<> urdist(0,1);
        NT M=1.0;
        while(true){
            r=urdist(var.rng);
            dis = (1.0-r)*l+r*u;
            r_val = M*urdist(var.rng);
            fn = std::exp(-a_i*dis*dis);
            if(r_val<fn){
                break;
            }
        }
    }

    return 1;
}


template <class T, class K>
int rand_gaussian_point_generator(T &P,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         NT a_i,
                         vars &var)  // constans for volume
{
    //std::cout<<"EDW2!"<<std::endl;
    int n = var.n;
    //bool birk = var.birk;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0,1);
    boost::random::uniform_int_distribution<> uidist(0,n-1);
    //std::uniform_real_distribution<NT> urdist = var.urdist;
    //std::uniform_int_distribution<int> uidist(0,n-1);

    std::vector<NT> lamdas(P.num_of_hyperplanes(),NT(0));
    //int rand_coord = rand()%n;
    //double kapa = double(rand())/double(RAND_MAX);
    int rand_coord = uidist(rng);
    double kapa = urdist(rng);
    Point p_prev = p;
    if(var.coordinate){
        //std::cout<<"[1a]P dim: "<<P.dimension()<<std::endl;
        gaussian_hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord,a_i,lamdas,var,true);
        //std::cout<<"[1b]P dim: "<<P.dimension()<<std::endl;
    }else
        gaussian_hit_and_run(p,P,a_i,var);

    for(int i=1; i<=rnum; ++i){

        for(int j=0; j<walk_len; ++j){
            int rand_coord_prev = rand_coord;
            //rand_coord = rand()%n;
            //kapa = double(rand())/double(RAND_MAX);
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            if(var.coordinate){
                //std::cout<<"[1c]P dim: "<<P.dimension()<<std::endl;
                gaussian_hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord_prev,a_i,lamdas,var,false);
            }else
                gaussian_hit_and_run(p,P,a_i,var);
        }
        randPoints.push_back(p);
        //if(birk) birk_sym(P,randPoints,p);
    }
    return 1;

    //if(rand_only) std::cout<<p<<std::endl;
    //if(print) std::cout<<"("<<i<<") Random point: "<<p<<std::endl;
}


template <class T>
int gaussian_hit_and_run(Point &p,
                T &P,
                NT a_i,
                vars &var)
{
    int n = var.n;
    double err = var.err;
    RNGType &rng = var.rng;
    //std::uniform_real_distribution<NT> &urdist = var.urdist;
    boost::random::uniform_real_distribution<> urdist(0,1);
    //std::uniform_real_distribution<NT> &urdist1 = var.urdist1;

    Point origin(n);

    Random_points_on_sphere_d<Point> gen (n, 1.0);
    Point l = gen.sample_point(rng);// - CGAL::Origin();
    //Point l2=origin;
    //Vector b1 = line_bisect(p,l,P,var,var2);
    //Vector b2 = line_bisect(p,-l,P,var,var2);
    //std::pair<Point,Point> ppair = P.line_intersect(p,l);
    std::pair<NT,NT> dbpair = P.line_intersect(p,l);
    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    NT dis;
    rand_exp_range(max_minus, min_plus, a_i, dis, var);
    p = (dis*l)+p;
    //Point b1 = (min_plus*l)+p;
    //Point b2 = (max_minus*l)+p;
    //Point b1 = ppair.first;// - origin;
    //Point b2 = ppair.second;// - origin;
    //std::cout<<"b1="<<b1<<"b2="<<b2<<std::endl;

    //NT lambda = urdist(rng);
    //p = (lambda*b1);
    //p=((1-lambda)*b2) + p;
    return 1;
}


//hit-and-run with orthogonal directions and update
template <class T>
int gaussian_hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             T &P,
                             int rand_coord,
                             int rand_coord_prev,
                             NT a_i,
                             std::vector<NT> &lamdas,
                             vars &var,
                             bool init)
{
    //std::cout<<"[1]P dim: "<<P.dimension()<<std::endl;
    std::pair<NT,NT> bpair;
    // EXPERIMENTAL
    //if(var.NN)
    //  bpair = P.query_dual(p,rand_coord);
    //else
    bpair = P.line_intersect_coord(p,p_prev,rand_coord,rand_coord_prev,lamdas,init);
    NT min_plus = bpair.first;
    NT max_minus = bpair.second;
    NT dis;
    rand_exp_range_coord(max_minus, min_plus, a_i, dis, var);
    p_prev = p;
    p.set_coord(rand_coord , p[rand_coord] + dis);
    return 1;
}


#endif