// VolEsti (volume computation and sampling library)

// Copyright (c) 2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Original C++ code from https://github.com/rzrsk/vaidya-walk by Raaz Dwivedi.

// Modified by Alexandros Manochis to be integrated in volesti, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

// The implemented random walk is presented in the paper of
// Y. Chen, R. Dwivedi, M. J. Wainwright and B. Yu,
// "Fast MCMC Sampling Algorithms on Polytopes",
// Journal of Machine Learning Research, 2018.

#ifndef PWALK_DIKIN_WALKER_HPP_
#define PWALK_DIKIN_WALKER_HPP_

#include <Eigen/Dense>
#include "math_functions.h"

/// @brief Class that defines the Dikin walk sampler
/// @tparam Dtype Number Type
template <typename Dtype>
class DikinWalker {
public:

    DikinWalker() {}

    DikinWalker(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1> &initialization, const Eigen::Matrix <Dtype,
    Eigen::Dynamic, Eigen::Dynamic> &cons_A, const Eigen::Matrix<Dtype, Eigen::Dynamic, 1> &cons_b, const Dtype r){
        nb_dim_ = cons_A.cols();
        nb_cons_ = cons_A.rows();
        nb_curr_samples_ = 1;
        initialization_ = initialization;
        cons_A_ = cons_A;
        cons_b_ = cons_b;
        curr_sample_ = initialization;
        r_ = r;
    }

    // getter for radius
    Dtype getRadius() {
        return r_;
    }

    // check whether a given point is in th polytope
    bool checkInPolytope(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1> &new_sample) {
        return (cons_A_ * new_sample - cons_b_).maxCoeff() < 0;
    }

    // getter for dimension
    int getNbDim() {
        return nb_dim_;
    }

    // getter for nb current samples
    int getNbCurrSamples() {
        return nb_curr_samples_;
    }

    // getter for current sample
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> &getCurrSample() {
        return curr_sample_;
    }

    void proposal(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
        Eigen::Matrix<Dtype, Eigen::Dynamic, 1> gaussian_step =
                        Eigen::Matrix<Dtype, Eigen::Dynamic, 1>::Zero(this->nb_dim_);
        sample_gaussian<Dtype>(this->nb_dim_, 0., 1., gaussian_step);

        // get hessian
        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess =
                        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
        sqrtInvHessBarrier(this->curr_sample_, new_sqrt_inv_hess);

        new_sample = this->curr_sample_ + r_ / std::sqrt(Dtype(this->nb_dim_))  * (new_sqrt_inv_hess * gaussian_step);
    }

    bool acceptRejectReverse(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
        // get hessian on x
        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess_x =
                        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
        sqrtInvHessBarrier(this->curr_sample_, new_sqrt_inv_hess_x);
        // get hessian on y
        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess_y =
                        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
        sqrtInvHessBarrier(new_sample, new_sqrt_inv_hess_y);

        Dtype scale = r_/std::sqrt(Dtype(this->nb_dim_));
        Dtype p_y_to_x = gaussian_density<Dtype>(this->curr_sample_, new_sample, new_sqrt_inv_hess_y.inverse()/scale);
        Dtype p_x_to_y = gaussian_density<Dtype>(new_sample, this->curr_sample_, new_sqrt_inv_hess_x.inverse()/scale);

        Dtype ar_ratio = std::min<Dtype>(1., p_y_to_x/p_x_to_y);

        Dtype random_num = rng_uniform<Dtype>(0., 1.);
        // lazy version of the walk
        if (random_num > ar_ratio) {
            return false;
        }

        return true;
    }

    bool doSample(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample, const Dtype lazy = Dtype(0.5)){
        proposal(new_sample);
        this->nb_curr_samples_ += 1;
        // for lazy markov chain
        Dtype random_num = rng_uniform<Dtype>(0., 1.);
        // check balance and check in polytope
        if (random_num < lazy && this->checkInPolytope(new_sample) && acceptRejectReverse(new_sample)){
            this->curr_sample_ = new_sample;
            return true;
        } else {
            new_sample = this->curr_sample_;
            return false;
        }
    }

    void sqrtInvHessBarrier(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample, Eigen::Matrix<Dtype,
                            Eigen::Dynamic, Eigen::Dynamic>& new_sqrt_inv_hess)
    {
        Eigen::Matrix<Dtype, Eigen::Dynamic, 1> inv_slack = (this->cons_b_ - this->cons_A_ * new_sample).cwiseInverse();

        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> half_hess = inv_slack.asDiagonal()* this->cons_A_;
        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_hess = half_hess.transpose() * half_hess;

        // compute eigenvectors and eigenvalues
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> > es(new_hess);

        Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> V = es.eigenvectors();
        Eigen::Matrix<Dtype, Eigen::Dynamic, 1> Dv = es.eigenvalues();
        new_sqrt_inv_hess = V * Dv.cwiseInverse().cwiseSqrt().asDiagonal() * V.transpose();
    }

private:
    Dtype r_;

    // Dimension
    int nb_dim_;
    // number of constraints
    int nb_cons_;
    // current sample size
    int nb_curr_samples_;
    // Initial vector
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> initialization_;
    // constraints matrix A
    Eigen::Matrix <Dtype, Eigen::Dynamic, Eigen::Dynamic> cons_A_;
    // constraint vector b
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> cons_b_;
    // Current vector
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> curr_sample_;
};


#endif // PWALK_DIKIN_WALKER_HPP_
