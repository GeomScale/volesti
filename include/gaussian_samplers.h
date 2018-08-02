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


template <typename FT>
FT eval_exp(Point p, FT a) {
    return std::exp(-a * p.squared_length());
}


template <typename FT>
FT get_max(Point l, Point u, FT a_i) {
    FT res;
    Point a = -1.0 * l;
    Point bef = u - l;
    Point b = (1.0 / std::sqrt((bef).squared_length())) * bef;
    Point z = (a.dot(b) * b) + l;
    FT low_bd = (l[0] - z[0]) / b[0];
    FT up_bd = (u[0] - z[0]) / b[0];
    if (low_bd * up_bd > 0) {
        //if(std::signbit(low_bd)==std::signbit(up_bd)){
        res = std::max(eval_exp(u, a_i), eval_exp(l, a_i));
    } else {
        res = eval_exp(z, a_i);
    }

    return res;
}


template <typename FT>
FT get_max_coord(FT l, FT u, FT a_i) {
    if (l < 0.0 && u > 0.0) {
        return 1.0;
    }
    return std::max(std::exp(-a_i * l * l), std::exp(-a_i * u * u));
}


template <typename FT>
void rand_exp_range(Point lower, Point upper, FT a_i, Point &p, vars_g &var) {
    FT r, r_val, fn;
    const FT tol = 0.00000001;
    Point bef = upper - lower;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    if (a_i > tol && std::sqrt(bef.squared_length()) >= (2.0 / std::sqrt(2.0 * a_i))) {
        boost::normal_distribution<> rdist(0, 1);
        Point a = -1.0 * lower;
        Point b = (1.0 / std::sqrt(bef.squared_length())) * bef;
        Point z = (a.dot(b) * b) + lower;
        FT low_bd = (lower[0] - z[0]) / b[0];
        FT up_bd = (upper[0] - z[0]) / b[0];
        while (true) {
            r = rdist(rng2);
            r = r / std::sqrt(2.0 * a_i);
            if (r >= low_bd && r <= up_bd) {
                break;
            }
        }
        p = (r * b) + z;

    } else {
        boost::random::uniform_real_distribution<> urdist(0, 1);
        FT M = get_max(lower, upper, a_i);
        while (true) {
            r = urdist(rng2);
            Point pef = r * upper;
            p = ((1.0 - r) * lower) + pef;
            r_val = M * urdist(var.rng);
            fn = eval_exp(p, a_i);
            if (r_val < fn) {
                break;
            }
        }
    }
}


template <typename FT>
FT rand_exp_range_coord(FT l, FT u, FT a_i, vars_g &var) {
    FT r, r_val, fn, dis;
    const FT tol = 0.00000001;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    if (a_i > tol && u - l >= 2.0 / std::sqrt(2.0 * a_i)) {
        boost::normal_distribution<> rdist(0, 1);
        while (true) {
            r = rdist(rng2);
            r = r / std::sqrt(2.0 * a_i);
            if (r >= l && r <= u) {
                break;
            }
        }
        dis = r;

    } else {
        boost::random::uniform_real_distribution<> urdist(0, 1);
        FT M = get_max_coord(l, u, a_i);
        while (true) {
            r = urdist(rng2);
            dis = (1.0 - r) * l + r * u;
            r_val = M * urdist(rng2);
            fn = std::exp(-a_i * dis * dis);
            if (r_val < fn) {
                break;
            }
        }
    }
    return dis;
}


template <class T, typename FT>
void gaussian_next_point(T &P,
                        Point &p,   // a point to start
                        Point &p_prev,
                        int &coord_prev,
                        int walk_len,
                        FT a_i,
                        std::vector<FT> &lamdas,
                        vars_g &var,
                        bool first) {
    int n = var.n, rand_coord;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    FT ball_rad = var.delta;

    if (var.coordinate && !var.ball_walk) {
        rand_coord = uidist(rng2);
        std::pair <FT, FT> bpair = P.line_intersect_coord(p,rand_coord, lamdas);
        //FT dis;
        FT dis = rand_exp_range_coord(p[rand_coord] + bpair.second, p[rand_coord] + bpair.first, a_i, var);
        p_prev = p;
        coord_prev = rand_coord;
        p.set_coord(rand_coord, dis);
        walk_len--;
    }

    for (int j = 0; j < walk_len; j++) {
        if (var.ball_walk) {
            gaussian_ball_walk(p, P, a_i, ball_rad, var);
        } else if (!var.coordinate) {
            gaussian_hit_and_run(p, P, a_i, var);
        } else {
            rand_coord = uidist(rng2);
            gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, a_i, lamdas, var);
            coord_prev = rand_coord;
        }
    }
}


template <class T, typename FT>
void gaussian_next_point(T &P,
                        Point &p,   // a point to start
                        Point &p_prev,
                        int &coord_prev,
                        int walk_len,
                        FT a_i,
                        std::vector<FT> &lamdas,
                        vars_g &var) {
    int n = var.n, rand_coord;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    FT ball_rad = var.delta;

    for (int j = 0; j < walk_len; j++) {
        if (var.ball_walk) {
            gaussian_ball_walk(p, P, a_i, ball_rad, var);
        } else if (!var.coordinate) {
            gaussian_hit_and_run(p, P, a_i, var);
        } else {
            rand_coord = uidist(rng2);
            gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, a_i, lamdas, var);
            coord_prev = rand_coord;
        }
    }
}


template <class T, class K, typename FT>
void rand_gaussian_point_generator(T &P,
                         Point &p,   // a point to start
                         int rnum,
                         int walk_len,
                         K &randPoints,
                         FT a_i,
                         vars_g &var)  // constans for volume
{
    int n = var.n;
    //bool birk = var.birk;
    //RNGType &rng = var.rng;
    RNGType &rng2 = var.rng;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    std::vector <FT> lamdas(P.num_of_hyperplanes(), FT(0));
    int rand_coord = uidist(rng2), coord_prev;
    FT ball_rad = var.delta;
    Point p_prev = p;

    if (var.coordinate && !var.ball_walk) {
        rand_coord = uidist(rng2);
        std::pair <FT, FT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        //FT dis;
        FT dis = rand_exp_range_coord(p[rand_coord] + bpair.second, p[rand_coord] + bpair.first, a_i, var);
        p_prev = p;
        coord_prev = rand_coord;
        p.set_coord(rand_coord, dis);
        for (int j = 0; j < walk_len - 1; ++j) {
            rand_coord = uidist(rng2);
            gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, a_i, lamdas, var);
            coord_prev = rand_coord;
        }
        randPoints.push_back(p);
        rnum--;
    }

    for (int i = 1; i <= rnum; ++i) {

        for (int j = 0; j < walk_len; ++j) {
            int rand_coord_prev = rand_coord;
            rand_coord = uidist(rng2);
            if (var.ball_walk) {
                gaussian_ball_walk(p, P, a_i, ball_rad, var);
            } else if (var.coordinate) {
                gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, a_i, lamdas, var);
            } else
                gaussian_hit_and_run(p, P, a_i, var);
        }
        randPoints.push_back(p);
    }
}


template <class T, typename FT>
void gaussian_hit_and_run(Point &p,
                T &P,
                FT a_i,
                vars_g &var) {
    int n = var.n;
    RNGType rng2 = var.rng;
    //Point l = get_point_on_Dsphere(n, 1.0);
    Point l = get_direction(n);
    std::pair <FT, FT> dbpair = P.line_intersect(p, l);

    FT min_plus = dbpair.first;
    FT max_minus = dbpair.second;
    Point upper = (min_plus * l) + p;
    Point lower = (max_minus * l) + p;

    rand_exp_range(lower, upper, a_i, p, var);
}


//hit-and-run with orthogonal directions and update
template <class T, typename FT>
void gaussian_hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             T &P,
                             int rand_coord,
                             int rand_coord_prev,
                             FT a_i,
                             std::vector<FT> &lamdas,
                             vars_g &var) {
    std::pair <FT, FT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    //FT dis;
    FT dis = rand_exp_range_coord(p[rand_coord] + bpair.second, p[rand_coord] + bpair.first, a_i, var);
    p_prev = p;
    p.set_coord(rand_coord, dis);
}


template <class T, typename FT>
void gaussian_ball_walk(Point &p,
              T &P,
              FT a_i,
              FT ball_rad,
              vars_g var) {
    int n = P.dimension();
    FT f_x, f_y, rnd;
    Point y = get_point_in_Dsphere(n, ball_rad);
    y = y + p;
    f_x = eval_exp(p, a_i);
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    boost::random::uniform_real_distribution<> urdist(0, 1);
    if (P.is_in(y) == -1) {
        f_y = eval_exp(y, a_i);
        rnd = urdist(rng2);
        if (rnd <= f_y / f_x) {
            p = y;
        }
    }
}

#endif
