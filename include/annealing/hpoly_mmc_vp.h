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

#ifndef HPOLY_MMC_VP_H
#define HPOLY_MMC_VP_H


template <class Polytope, class Ball, class Parameters>
void enclosing_ball(Polytope &P, Ball &B0, Parameters &var) {

    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Polytope::NT NT;
    typedef typename Polytope::PolytopePoint Point;
    unsigned int n = P.dimension();
    std::vector<NT> vec(n,0.0);

    MT V = P.get_mat();
    //P.print();
    int k = V.rows();
    Point temp;
    NT rad = 0.0, nr;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            vec[j] = V(i,j);
        }
        temp = Point(n, vec.begin(), vec.end());
        nr = std::sqrt(temp.squared_length());
        if ( nr > rad) {
            rad = nr;
        }
    }
    rad = rad*1.1;
    //std::cout<<"radius of minim ball = "<<rad<<std::endl;
    //Point xc = Point(n);
    //std::vector<Ball> S0;
    B0 = Ball(Point(n), rad*rad);
    //ConvSet.push_back(Interballs(P.dimension(), vecBalls));

}


template <class HPolyBall, class VPolytope, class Hpolytope, class Point, class Parameters>
void add_new_facet(HPolyBall &PB, VPolytope &VP, Hpolytope &HP, Point &q, Parameters &var){

    unsigned int n = var.n, count = 0;
    typedef typename VPolytope::NT NT;
    typedef typename VPolytope::VT VT;
    typedef typename VPolytope::MT MT;
    Point center(n);
    std::vector <NT> lambdas1(VP.num_of_vertices()), lambdas2(VP.num_of_vertices());
    NT z0 = 1.0;
    //std::cout<<P.is_in(center)<<std::endl;
    Point v = q * (1.0 / std::sqrt(q.squared_length()));
    VP.intersect_double_line_Vpoly_return(center, v, lambdas1, lambdas2);
    MT Mat(n,n), V = VP.get_mat();
    VT p(n);

    for (int j = 0; j < VP.num_of_vertices(); ++j) {
        if (lambdas1[j] > 0.0) {
            Mat.row(count) = V.row(j);
            count++;
        } else {
            p = V.row(j);
        }
    }
    //std::cout<<"count = "<<count<<"\n Mat = "<<Mat<<"\n"<<std::endl;
    VT a = Mat.colPivHouseholderQr().solve(VT::Ones(n));
    if (a.dot(p) > 1.0) {
        a = -a;
        z0 = -1.0;
    }
    PB.add_facet(a, z0);
    HP.add_facet(a, z0);
    //HP.print();

}

template <class Vpolytope, class Hpolytope, class PolyBall, class Parameters>
void construct_hpoly(Vpolytope &VP, Hpolytope &HP, PolyBall &PB, int walkL, int k, Parameters &var){

    unsigned int n = var.n, coord_prev;
    typedef typename Hpolytope::PolytopePoint Point;
    typedef typename Vpolytope::NT NT;
    typedef typename Vpolytope::VT VT;
    typedef typename Vpolytope::MT MT;

    Point p(n), p_prev(n);
    //NT lambda;
    //std::vector <NT> lamdas(HP.num_of_hyperplanes(), NT(0)), Av(HP.num_of_hyperplanes(), NT(0));

    //uniform_first_point(PB, p, p_prev, coord_prev, walkL, lamdas, Av, lambda, var);

    for (int i = 0; i < k; ++i) {

        //std::cout<<"sampling new point, walkL ="<<walkL<<std::endl;
        p = Point(n);
        do {
            for (int i = 0; i < walkL; ++i) {
                hit_and_run(p, PB, var);
            }
        }while(VP.is_in(p)==-1);
        //std::cout<<"p in VP = "<<VP.is_in(p)<<std::endl;
        //uniform_next_point(PB, p, p_prev, coord_prev, walkL, lamdas, Av, lambda, var);
        //std::cout<<"adding new facet"<<std::endl;
        add_new_facet(PB, VP, HP, p, var);

    }

}


#endif
