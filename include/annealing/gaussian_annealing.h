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
    typedef typename std::vector<FT>::iterator viterator;

    int i=0;
    viterator vecit = vec.begin();
    for( ; vecit!=vec.end(); vecit++, i++){
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }

    return std::pair<FT,FT> (mean, variance);
}


// Compute the first variance a_0 for the starting gaussian
template <class T1, typename FT>
void get_first_gaussian(T1 &K, FT radius, FT frac, const vars_g var, FT &error, std::vector<FT> &a_vals) {

    int m = K.num_of_hyperplanes(), dim = var.n, i=0;
    const int maxiter = 10000;
    const FT tol = 0.0000001;
    FT sum, lower = 0.0, upper = 1.0, mid;
    std::vector <FT> dists(m, 0);
    Eigen::MatrixXd A = K.get_eigen_mat();
    Eigen::VectorXd b = K.get_eigen_vec();
    typedef typename std::vector<FT>::iterator viterator;

    viterator disit = dists.begin();
    for ( ; disit!=dists.end(); disit++, i++)
        *disit = b(i)/A.row(i).norm();

    // Compute an upper bound for a_0
    for (i= 1; i <= maxiter; ++i) {
        sum = 0.0;
        for (viterator it = dists.begin(); it != dists.end(); ++it) {
            sum += std::exp(-upper * std::pow(*it, 2.0)) / (2.0 * (*it) * std::sqrt(M_PI * upper));
        }

        if (sum > frac * error) {
            upper = upper * 10;
        } else {
            break;
        }
    }

    if (i == maxiter) {
        std::cout << "Cannot obtain sharp enough starting Gaussian" << std::endl;
        exit(-1);
    }

    //get a_0 with binary search
    while (upper - lower > tol) {
        mid = (upper + lower) / 2.0;
        sum = 0.0;
        for (viterator it = dists.begin(); it != dists.end(); ++it) {
            sum += std::exp(-mid * std::pow(*it, 2)) / (2 * (*it) * std::sqrt(M_PI * mid));
        }

        if (sum < frac * error) {
            upper = mid;
        } else {
            lower = mid;
        }
    }

    a_vals.push_back((upper + lower) / 2.0);
    error = (1.0 - frac) * error;
}


// Compute a_{i+1} when a_i is given
template <class T1, typename FT>
FT get_next_gaussian(T1 K,Point &p, FT a, int N, FT ratio, FT C, vars_g var){

    FT last_a = a, last_ratio = 0.1;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    FT k = 1.0;
    const FT tol = 0.00001;
    bool done=false, print=var.verbose;
    std::vector<FT> fn(N,0);
    std::list<Point> randPoints;
    typedef typename std::vector<FT>::iterator viterator;

    //sample N points using hit and run or ball walk
    rand_gaussian_point_generator(K, p, N, var.walk_steps, randPoints, last_a, var);

    viterator fnit;
    while(!done){
        a = last_a*std::pow(ratio,k);

        fnit = fn.begin();
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, fnit++){
            *fnit = eval_exp(*pit,a)/eval_exp(*pit, last_a);
        }
        std::pair<FT,FT> mv = getMeanVariance(fn);

        // Compute a_{i+1}
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
    return last_a*std::pow(ratio, k);
}


// Compute the sequence of spherical gaussians
template <class T1, typename FT>
void get_annealing_schedule(T1 K, FT radius, FT ratio, FT C, FT frac, int N, vars_g var, FT &error, std::vector<FT> &a_vals){
    bool print=var.verbose;
    // Compute the first gaussian
    get_first_gaussian(K, radius, frac, var, error, a_vals);
    if(print) std::cout<<"first gaussian computed\n"<<std::endl;

    FT a_stop = 0.0, curr_fn = 2.0, curr_its = 1.0, next_a;
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
    while (true) {

        if (var.ball_walk) {
            if (var.delta < 0.0) {
                var.delta = 4.0 * var.che_rad / std::sqrt(std::max(1.0, a_vals[it]) * FT(n));
            }
        }
        // Compute the next gaussian
        next_a = get_next_gaussian(K, p, a_vals[it], N, ratio, C, var);

        curr_fn = 0;
        curr_its = 0;
        std::fill(lamdas.begin(), lamdas.end(), FT(0));
        steps = totalSteps;

        if (var.coordinate && !var.ball_walk){
            gaussian_first_coord_point(K, p, p_prev, coord_prev, var.walk_steps, a_vals[it], lamdas, var);
            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
            steps--;
        }


        // Compute some ratios to decide if this is the last gaussian
        for (int j = 0; j < steps; j++) {
            gaussian_next_point(K, p, p_prev, coord_prev, var.walk_steps, a_vals[it], lamdas, var);
            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
        }

        // Remove the last gaussian.
        // Set the last a_i equal to zero
        if (next_a>0 && curr_fn/curr_its>(1.0+tol)) {
            a_vals.push_back(next_a);
            it ++;
        } else if (next_a <= 0) {
            a_vals.push_back(a_stop);
            it++;
            break;
        } else {
            a_vals[it] = a_stop;
            break;
        }
    }

}


#endif
