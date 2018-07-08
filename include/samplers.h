// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

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

#ifndef RANDOM_SAMPLERS_H
#define RANDOM_SAMPLERS_H


template <class P>
class Random_points_on_sphere_d
{
public:
    typedef typename P::FT 	FT;
    int d;
    FT r;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(std::chrono::system_clock::now().time_since_epoch().count());
    //rn generator;
    //std::default_random_engine generator;   
    //std::normal_distribution<FT> distribution = std::normal_distribution<FT>(FT(0),FT(1));
    boost::normal_distribution<> distribution = boost::normal_distribution<>(0,1);
    //typedef std::vector<K> coeff;
    //coeff coeffs;
    
    Random_points_on_sphere_d() {}
    
    Random_points_on_sphere_d(int dim, FT radius){
        d=dim;
        r=radius;
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


// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //

template <class T, class K>
int rand_point_generator(T &P,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         vars &var  // constans for volume
                         )
{
	//std::cout<<"EDW2!"<<std::endl;
    int n = var.n;
    bool birk = var.birk;
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
        hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord,kapa,lamdas,var,var,true);
        //std::cout<<"[1b]P dim: "<<P.dimension()<<std::endl;
    }else
        hit_and_run(p,P,var,var);

    for(int i=1; i<=rnum; ++i){

        for(int j=0; j<walk_len; ++j){
            int rand_coord_prev = rand_coord;
            //rand_coord = rand()%n;
            //kapa = double(rand())/double(RAND_MAX);
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            if(var.coordinate){
				//std::cout<<"[1c]P dim: "<<P.dimension()<<std::endl;
                hit_and_run_coord_update(p,p_prev,P,rand_coord,rand_coord_prev,kapa,lamdas,var,var,false);
            }else
                hit_and_run(p,P,var,var);
        }
        randPoints.push_back(p);
        if(birk) birk_sym(P,randPoints,p);
    }

    //if(rand_only) std::cout<<p<<std::endl;
    //if(print) std::cout<<"("<<i<<") Random point: "<<p<<std::endl;
}



template <class T, class K>
int rand_point_generator(BallIntersectPolytope<T> &PBLarge,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         BallIntersectPolytope<T> &PBSmall,
                         int &nump_PBSmall,
                         vars &var  // constans for volume
                         )
{
	//std::cout<<"EDW!: "<<std::endl;
    int n = var.n;
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist(0,1);
    boost::random::uniform_int_distribution<> uidist(0,n-1);
    //std::uniform_real_distribution<NT> urdist = var.urdist;
    //std::uniform_int_distribution<int> uidist(0,n-1);

    std::vector<NT> lamdas(PBLarge.num_of_hyperplanes(),NT(0));
    //int rand_coord = rand()%n;
    //double kapa = double(rand())/double(RAND_MAX);
    int rand_coord = uidist(rng);
    double kapa = urdist(rng);
    Point p_prev = p;

    if(var.coordinate)
        hit_and_run_coord_update(p,p_prev,PBLarge,rand_coord,rand_coord,kapa,lamdas,var,var,true);
    else
        hit_and_run(p,PBLarge,var,var);

    for(int i=1; i<=rnum; ++i){
        for(int j=0; j<walk_len; ++j){
            int rand_coord_prev = rand_coord;
            //rand_coord = rand()%n;
            //kapa = double(rand())/double(RAND_MAX);
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            if(var.coordinate)
                hit_and_run_coord_update(p,p_prev,PBLarge,rand_coord,rand_coord_prev,kapa,lamdas,var,var,false);
            else
                hit_and_run(p,PBLarge,var,var);
        }
        if(PBSmall.second().is_in(p) == -1){//is in
            randPoints.push_back(p);
            ++nump_PBSmall;
        }
    }
    //if(rand_only) std::cout<<p<<std::endl;
    //if(print) std::cout<<"("<<i<<") Random point: "<<p<<std::endl;
}

// ----- HIT AND RUN FUNCTIONS ------------ //

//hit-and-run with random directions
template <class T>
int hit_and_run(Point &p,
                T &P,
                vars &var,
                vars &var2)
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
    std::pair<Point,Point> ppair = P.line_intersect(p,l);
    Point b1 = ppair.first;// - origin;
    Point b2 = ppair.second;// - origin;
    //std::cout<<"b1="<<b1<<"b2="<<b2<<std::endl;
    NT lambda = urdist(rng);
    p = (lambda*b1);
    p=((1-lambda)*b2) + p;
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
                             vars &var2,
                             bool init)
{	
	//std::cout<<"[1]P dim: "<<P.dimension()<<std::endl;
    std::pair<NT,NT> bpair;
    // EXPERIMENTAL
    //if(var.NN)
    //  bpair = P.query_dual(p,rand_coord);
    //else
    bpair = P.line_intersect_coord(p,p_prev,rand_coord,rand_coord_prev,lamdas,init);
    //std::cout<<"[2]P dim: "<<P.dimension()<<std::endl;
    //std::cout<<"original:"<<bpair.first<<" "<<bpair.second<<std::endl;
    //std::cout<<"-----------"<<std::endl;
    //TODO: only change one coordinate of *r* avoid addition + construction
    //std::vector<NT> v(P.dimension(),NT(0));
    //v[rand_coord] = bpair.first + kapa * (bpair.second - bpair.first);
    //Point vp(P.dimension(),v.begin(),v.end());
    p_prev = p;
    //std::cout<<"v dim: "<<v.size()<<std::endl;
   // std::cout<<"P dim: "<<P.dimension()<<std::endl;
    //p = p + vp;
    p.set_coord(rand_coord , p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
    return 1;
}



#endif //RANDOM_SAMPLERS_H
