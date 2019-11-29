// VolEsti (volume computation and sampling library)

// Copyright (c) 2019 Vissarion Fisikopoulos, Apostolos Chalkis

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

#ifndef HPOLY_MMC_PROJ_H
#define HPOLY_MMC_PROJ_H


template <class HPolyBall, class VPolytope, class Hpolytope, class Point, class Parameters>
void add_new_facet2(HPolyBall &PB, VPolytope &VP, Hpolytope &HP, Point &q, Parameters &var){

    unsigned int _d = var.n;
    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::MT MT;
    typedef boost::mt19937 RNGType;

    Point center(_d);
    int count = 0, k =VP.get_T().cols();
    std::vector <NT> lambdas1(k), lamdas2(k);//, lambdas2(VP.num_of_vertices());
    NT z0 = 1.0;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    //v = -1.0*v;
    VP.intersect_double_line_Vpoly_return(center, v, lambdas1, lamdas2);

    NT sum, e = 0.0000000001;
    MT A = VP.get_mat();
    VT b = VP.get_vec();
    MT T = VP.get_T();
    MT Fmat(k - _d +1, k), HypMat(_d, _d);
    VT rand_point(k), beq(k-_d+1);

    for (int i = 0; i < A.rows(); ++i) {
        sum = 0.0;
        for (int j = 0; j < k; ++j) {
            sum += A(i,j)* lambdas1[j];
        }
        if (std::abs(b(i) - sum)  < e*std::abs(b(i)) && std::abs(b(i) - sum)  < e*std::abs(sum)) {
            //std::cout<<"a_"<<i<<"x = "<<sum<<" b("<<i<<") = "<<b(i)<<std::endl;
            Fmat.row(count) = A.row(i);
            beq(count) = b(i);
            count++;
        }

    }
    //std::cout<<"rows equal to b = "<<count<<std::endl;
    //std::cout<<"\n"<<std::endl;

    Eigen::FullPivLU<MT> lu(Fmat);
    MT NS = lu.kernel();

    VT x_pr = Fmat.colPivHouseholderQr().solve(beq);
    VT x0 = T * x_pr;
    MT TT = T * NS;
    HypMat.row(0) = x0.transpose();

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    NT normal;
    for (int l = 0; l < _d-1; ++l) {
        normal = 0.0;
        for (unsigned int i = 0; i < k; i++) {
            rand_point(i) = rdist(rng);
            normal += rand_point(i) * rand_point(i);
        }
        normal = 1.0 / std::sqrt(normal);
        rand_point = rand_point * normal;
        HypMat.row(l+1) = (x0 + TT * rand_point).transpose();
        //std::cout<<"Hypmat = "<<HypMat<<"\n(x0 + TT * rand_point).transpose() ="<<(x0 + TT * rand_point).transpose()<<std::endl;
    }

    VT a = HypMat.colPivHouseholderQr().solve(VT::Ones(_d));
    //z0 = z0/a.norm();
    //a = a/a.norm();
    //std::cout<<"HypMat = \n"<<HypMat<<"\n"<<std::endl;
    //std::cout<<"a = "<<a(0)<<std::endl;
    //sum = 0.0;
    //for (int i = 0; i < _d; ++i) sum += a(i)*p[i];
    //if(sum<0.0){
        //z0 = -1.0;
        //a = -1.0*a;
    //}


    PB.add_facet(a, z0);
    HP.add_facet(a, z0);
    //HP.print();

}

template <class Vpolytope, class Hpolytope, class PolyBall, class Parameters>
void construct_hpoly2(Vpolytope &VP, Hpolytope &HP, PolyBall &PB, int walkL, int k, Parameters &var){

    unsigned int n = var.n, coord_prev;
    typedef typename Hpolytope::PolytopePoint Point;
    typedef typename Vpolytope::NT NT;
    typedef typename Vpolytope::VT VT;
    typedef typename Vpolytope::MT MT;

    Point p(n), p_prev(n);
    //NT lambda;
    //std::vector <NT> lamdas(HP.num_of_hyperplanes(), NT(0)), Av(HP.num_of_hyperplanes(), NT(0));

    //uniform_first_point(PB, p, p_prev, coord_prev, walkL, lamdas, Av, lambda, var);
    std::list<Point> randPoints;
    int counter=1;
    NT rad =-1.0;
    std::pair<Point , NT> pair_rad;
    while (counter<=k){

        //std::cout<<"sampling new point, walkL ="<<walkL<<std::endl;
        p = Point(n);
        do {
            rand_point_generator(PB, p, 1, walkL, randPoints, var);
            //std::cout<<"p in VP = "<<VP.is_in(p)<<std::endl;
        //for (int i = 0; i < walkL; ++i) {
        //hit_and_run(p, PB, var);
        //}
        }while(VP.is_in(p)==-1);
        //std::cout<<"p in VP = "<<VP.is_in(p)<<std::endl;
        //uniform_next_point(PB, p, p_prev, coord_prev, walkL, lamdas, Av, lambda, var);
        //std::cout<<"adding new facet"<<std::endl;
        add_new_facet2(PB, VP, HP, p, var);
        randPoints.clear();
        counter++;
        if (counter>k) {
            pair_rad = HP.ComputeInnerBall();
            rad = pair_rad.second;
        }

    }
    std::cout<<"\nChebychev radius in hpoly construction = "<<rad<<", num_of_hyps = "<<HP.num_of_hyperplanes()<<"\n"<<std::endl;

}


#endif
