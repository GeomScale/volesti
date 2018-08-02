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


#ifndef GAUSSIAN_ANNEALING_H
#define GAUSSIAN_ANNEALING_H

#include <complex>


//An implementation of Welford's algorithm for mean and variance.
template <typename FT>
std::pair<FT,FT> getMeanVariance(std::vector<FT>& vec) {
    FT mean = 0, M2 = 0, variance = 0, delta;

    int i=0;
    typename std::vector<FT>::iterator vecit = vec.begin();
    for( ; vecit!=vec.end(); vecit++, i++){
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }

    return std::pair<FT,FT> (mean, variance);
}



template <class T1, typename FT>
void get_first_gaussian(T1 &K, FT radius, FT frac, const vars_g var, FT &error, std::vector<FT> &a_vals) {

    int m = K.num_of_hyperplanes(), dim = var.n;
    unsigned int iterations = 0;
    const int maxiter = 10000;
    const FT tol = 0.0000001;
    FT sum, lower = 0.0, upper = 1.0, sigma_sqd, t, mid;
    std::vector <FT> dists(m, 0);
    for (int i = 0; i < m; i++) {
        sum = 0.0;
        for (int j = 1; j < dim + 1; j++) {
            sum += K.get_coeff(i, j) * K.get_coeff(i, j);
        }
        dists[i] = K.get_coeff(i, 0) / std::sqrt(sum);
    }

    while (iterations < maxiter) {
        iterations += 1;
        sum = 0.0;
        for (typename std::vector<FT>::iterator it = dists.begin(); it != dists.end(); ++it) {
            sum += std::exp(-upper * std::pow(*it, 2.0)) / (2.0 * (*it) * std::sqrt(M_PI * upper));
        }

        sigma_sqd = 1 / (2.0 * upper);

        if (sum > frac * error) {
            upper = upper * 10;
        } else {
            break;
        }
    }

    if (iterations == maxiter) {
        std::cout << "Cannot obtain sharp enough starting Gaussian" << std::endl;
        exit(-1);
    }

    //get a_0 with binary search
    while (upper - lower > tol) {
        mid = (upper + lower) / 2.0;
        sum = 0.0;
        for (typename std::vector<FT>::iterator it = dists.begin(); it != dists.end(); ++it) {
            sum += std::exp(-mid * std::pow(*it, 2)) / (2 * (*it) * std::sqrt(M_PI * mid));
        }

        sigma_sqd = 1.0 / (2.0 * mid);

        if (sum < frac * error) {
            upper = mid;
        } else {
            lower = mid;
        }
    }

    a_vals.push_back((upper + lower) / 2.0);
    error = (1.0 - frac) * error;
}


template <class T1, typename FT>
void get_next_gaussian(T1 K,Point &p, FT a, int N, FT ratio, FT C, vars_g var, std::vector<FT> &a_vals){

    FT last_a = a, last_ratio = 0.1;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    FT k = 1.0;
    const FT tol = 0.00001;
    bool done=false, print=var.verbose;
    std::vector<FT> fn(N,0);

    //sample N points using hit and run
    std::list<Point> randPoints;
    if (var.ball_walk) {
        if (var.delta < 0.0) {
            var.delta = 4.0 * var.che_rad / std::sqrt(std::max(1.0, last_a) * FT(var.n));
        }
    }
    rand_gaussian_point_generator(K, p, N, var.walk_steps, randPoints, last_a, var);
    typename std::vector<FT>::iterator fnit;

    while(!done){
        a = last_a*std::pow(ratio,(FT(k)));

        fnit = fn.begin();
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, fnit++){
            *fnit = eval_exp(*pit,a)/eval_exp(*pit, last_a);
        }
        std::pair<FT,FT> mv = getMeanVariance(fn);

        if(mv.second/(mv.first * mv.first)>=C || mv.first/last_ratio<1.0+tol){
            if(k!=1.0){
                k=k/2;
            }
            done=true;
        }else{
            k=2*k;
        }
        last_ratio = mv.first;
    }
    a_vals.push_back(last_a*std::pow(ratio, k ) );
}


template <class T1, typename FT>
void get_annealing_schedule(T1 K, FT radius, FT ratio, FT C, FT frac, int N, vars_g var, FT &error, std::vector<FT> &a_vals){
    bool print=var.verbose;
    get_first_gaussian(K, radius, frac, var, error, a_vals);
    if(print) std::cout<<"first gaussian computed\n"<<std::endl;
    FT a_stop = 0.0, curr_fn = 2.0, curr_its = 1.0;
    const FT tol = 0.001;
    int it = 0, n = var.n, steps, coord_prev;
    const int totalSteps= ((int)150/error)+1;
    std::list<Point> randPoints;

    if(a_vals[0]<a_stop) {
        a_vals[0] = a_stop;
    }

    Point p(n);

    if(print) std::cout<<"Computing the sequence of gaussians..\n"<<std::endl;

    Point p_prev=p;
    std::vector<FT> lamdas(K.num_of_hyperplanes(),NT(0));
    while (curr_fn/curr_its>(1.0+tol) && a_vals[it]>=a_stop) {
        get_next_gaussian(K, p, a_vals[it], N, ratio, C, var, a_vals);
        it++;

        curr_fn = 0;
        curr_its = 0;
        std::fill(lamdas.begin(), lamdas.end(), FT(0));
        steps = totalSteps;

        if (var.coordinate && !var.ball_walk){
            gaussian_next_point(K, p, p_prev, coord_prev, var.walk_steps, a_vals[it - 1], lamdas, var, first_coord_point);
            curr_its += 1.0;
            curr_fn += eval_exp(p, a_vals[it]) / eval_exp(p, a_vals[it - 1]);
            steps--;
        }
        if (var.ball_walk) {
            if (var.delta < 0.0) {
                var.delta = 4.0 * var.che_rad / std::sqrt(std::max(1.0, a_vals[it - 1]) * FT(n));
            }
        }

        for (int j = 0; j < steps; j++) {
            gaussian_next_point(K, p, p_prev, coord_prev, var.walk_steps, a_vals[it - 1], lamdas, var);
            curr_its += 1.0;
            curr_fn += eval_exp(p, a_vals[it]) / eval_exp(p, a_vals[it - 1]);
        }
    }
    if (a_vals[it]>a_stop) {
        a_vals.pop_back();
        a_vals[it - 1] = a_stop;

    }else {
        a_vals[it] = a_stop;
    }
}


#endif
