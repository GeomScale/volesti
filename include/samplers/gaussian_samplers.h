// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

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


// evaluate the pdf of point p
template <typename Point, typename NT>
NT eval_exp(Point &p, const NT &a) {
    return std::exp(-a * p.squared_length());
}


template <typename Point, typename NT>
NT get_max(Point &l, Point &u, const NT &a_i) {
    NT res;
    Point a = -1.0 * l, bef = u - l;
    Point b = (1.0 / std::sqrt((bef).squared_length())) * bef;
    Point z = (a.dot(b) * b) + l;
    NT low_bd = (l[0] - z[0]) / b[0], up_bd = (u[0] - z[0]) / b[0];
    if (low_bd * up_bd > 0) {
        //if(std::signbit(low_bd)==std::signbit(up_bd)){
        res = std::max(eval_exp(u, a_i), eval_exp(l, a_i));
    } else {
        res = eval_exp(z, a_i);
    }

    return res;
}


template <typename NT>
NT get_max_coord(const NT &l, const NT &u, const NT &a_i) {
    if (l < 0.0 && u > 0.0) {
        return 1.0;
    }
    return std::max(std::exp(-a_i * l * l), std::exp(-a_i * u * u));
}


// Pick a point from the distribution exp(-a_i||x||^2) on the chord
template <typename Parameters, typename Point, typename NT>
void rand_exp_range(Point &lower, Point & upper, const NT &a_i, Point &p, Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    NT r, r_val, fn;
    const NT tol = 0.00000001;
    Point bef = upper - lower;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
    if (a_i > tol && std::sqrt(bef.squared_length()) >= (2.0 / std::sqrt(2.0 * a_i))) {
        boost::normal_distribution<> rdist(0, 1);
        Point a = -1.0 * lower;
        Point b = (1.0 / std::sqrt(bef.squared_length())) * bef;
        Point z = (a.dot(b) * b) + lower;
        NT low_bd = (lower[0] - z[0]) / b[0];
        NT up_bd = (upper[0] - z[0]) / b[0];
        while (true) {
            r = rdist(rng2);
            r = r / std::sqrt(2.0 * a_i);
            if (r >= low_bd && r <= up_bd) {
                break;
            }
        }
        p = (r * b) + z;

    // select using rejection sampling from a bounding rectangle
    } else {
        boost::random::uniform_real_distribution<> urdist(0, 1);
        NT M = get_max(lower, upper, a_i);
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


// Pick a point from the distribution exp(-a_i||x||^2) on the coordinate chord
template <typename Parameters, typename NT>
NT rand_exp_range_coord(const NT &l, const NT &u, const NT &a_i, Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    NT r, r_val, fn, dis;
    const NT tol = 0.00000001;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
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

    // select using rejection sampling from a bounding rectangle
    } else {
        boost::random::uniform_real_distribution<> urdist(0, 1);
        NT M = get_max_coord(l, u, a_i);
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


// compute the first coordinate point
template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void gaussian_first_coord_point(Polytope &P,
                         Point &p,   // a point to start
                         Point &p_prev, // previous point
                         unsigned int &coord_prev, // previous coordinate ray
                         unsigned int walk_len, // number of steps for the random walk
                         const NT &a_i,
                         VT &lamdas,
                         Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n, rand_coord;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;

    rand_coord = uidist(rng2);
    std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
    NT dis = rand_exp_range_coord(p[rand_coord] + bpair.second, p[rand_coord] + bpair.first, a_i, var);
    p_prev = p;
    coord_prev = rand_coord;
    p.set_coord(rand_coord, dis);
    walk_len--;

    for (unsigned int j = 0; j < walk_len; j++) {
        rand_coord = uidist(rng2);
        gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, a_i, lamdas, var);
        coord_prev = rand_coord;
    }
}


// Compute the next point with target distribution the gaussian
template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void gaussian_next_point(Polytope &P,
                        Point &p,   // a point to start
                        Point &p_prev, // previous point
                        unsigned int &coord_prev, // previous coordinate ray
                        const unsigned int walk_len, // number of steps for the random walk
                        const NT &a_i,
                        VT &lamdas,
                        Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n, rand_coord;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    NT ball_rad = var.delta;

    for (unsigned int j = 0; j < walk_len; j++) {
        if (var.ball_walk) {
            gaussian_ball_walk(p, P, a_i, ball_rad, var);
        } else if (var.cdhr_walk) {
            rand_coord = uidist(rng2);
            gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, a_i, lamdas, var);
            coord_prev = rand_coord;
        } else {
            gaussian_hit_and_run(p, P, a_i, var);
        }
    }
}


// Sample N points with target distribution the gaussian
template <typename Polytope, typename Parameters, typename Point, typename PointList, typename NT>
void rand_gaussian_point_generator(Polytope &P,
                         Point &p,   // a point to start
                         unsigned int rnum,   // number of points to sample
                         const unsigned int walk_len,  // number of stpes for the random walk
                         PointList &randPoints,  // list to store the sampled points
                         const NT &a_i,
                         Parameters const& var)  // constans for volume
{
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    RNGType &rng2 = var.rng;
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    VT lamdas;
    lamdas.setZero(P.num_of_hyperplanes());

    unsigned int rand_coord = uidist(rng2), coord_prev, rand_coord_prev;
    NT ball_rad = var.delta;
    Point p_prev = p;

    if (var.cdhr_walk) {
        rand_coord = uidist(rng2);
        std::pair <NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        NT dis = rand_exp_range_coord(p[rand_coord] + bpair.second, p[rand_coord] + bpair.first, a_i, var);
        p_prev = p;
        coord_prev = rand_coord;
        p.set_coord(rand_coord, dis);
        for (unsigned int j = 0; j < walk_len - 1; ++j) {
            rand_coord = uidist(rng2);
            gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, coord_prev, a_i, lamdas, var);
            coord_prev = rand_coord;
        }
        randPoints.push_back(p);
        rnum--;
    }

    for (unsigned  int i = 1; i <= rnum; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            if (var.ball_walk) {
                gaussian_ball_walk(p, P, a_i, ball_rad, var);
            } else if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng2);
                gaussian_hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, a_i, lamdas, var);
            } else
                gaussian_hit_and_run(p, P, a_i, var);
        }
        randPoints.push_back(p);
    }
}


// hit-and-run with random directions and update
template <typename Polytope, typename Parameters, typename Point, typename NT>
void gaussian_hit_and_run(Point &p,
                Polytope &P,
                const NT &a_i,
                Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = var.n;
    RNGType rng2 = var.rng;
    Point l = get_direction<RNGType, Point, NT>(n);
    std::pair <NT, NT> dbpair = P.line_intersect(p, l);

    NT min_plus = dbpair.first;
    NT max_minus = dbpair.second;
    Point upper = (min_plus * l) + p;
    Point lower = (max_minus * l) + p;

    rand_exp_range(lower, upper, a_i, p, var);
}


// hit-and-run with orthogonal directions and update
template <typename Polytope, typename Parameters, typename Point, typename NT, typename VT>
void gaussian_hit_and_run_coord_update(Point &p,
                             Point &p_prev,
                             Polytope &P,
                             unsigned int rand_coord,
                             unsigned int rand_coord_prev,
                             const NT &a_i,
                             VT &lamdas,
                             Parameters const& var) {
    std::pair <NT, NT> bpair = P.line_intersect_coord(p, p_prev, rand_coord, rand_coord_prev, lamdas);
    NT dis = rand_exp_range_coord(p[rand_coord] + bpair.second, p[rand_coord] + bpair.first, a_i, var);
    p_prev = p;
    p.set_coord(rand_coord, dis);
}



// ball walk and update
template <typename Polytope, typename Parameters, typename Point, typename NT>
void gaussian_ball_walk(Point & p,
                        Polytope & P,
                        NT const& a_i,
                        NT const& ball_rad,
                        Parameters const& var) {
    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    NT f_x, f_y, rnd;
    Point y = get_point_in_Dsphere<RNGType, Point>(n, ball_rad);
    y += p;
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
