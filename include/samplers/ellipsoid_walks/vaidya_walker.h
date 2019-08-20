// Code from https://github.com/rzrsk/vaidya-walk

// Modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

#ifndef PWALK_VAIDYA_WALKER_HPP_
#define PWALK_VAIDYA_WALKER_HPP_

#include <Eigen/Dense>
#include "math_functions.h"
#include "walker.h"

//namespace pwalk {

template <typename Dtype>
class VaidyaWalker: public Walker<Dtype> {
public:
  VaidyaWalker(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& initialization, const Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& cons_A, const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& cons_b, const Dtype r) : Walker<Dtype>(initialization, cons_A, cons_b), r_(r){}

  // getter for radius
  Dtype getRadius() {
    return r_;
  }

  void proposal(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> gaussian_step = Eigen::Matrix<Dtype, Eigen::Dynamic, 1>::Zero(this->nb_dim_);
    sample_gaussian<Dtype>(this->nb_dim_, 0., 1., gaussian_step);

    // get hessian
    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess = Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
    sqrtInvHessBarrier(this->curr_sample_, new_sqrt_inv_hess);

    new_sample = this->curr_sample_ + r_ / std::sqrt(std::sqrt(Dtype(this->nb_dim_)*Dtype(this->nb_cons_)))  * (new_sqrt_inv_hess * gaussian_step);
  }

  bool acceptRejectReverse(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
    // get hessian on x
    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess_x = Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
    sqrtInvHessBarrier(this->curr_sample_, new_sqrt_inv_hess_x);
    // get hessian on y
    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_sqrt_inv_hess_y = Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->nb_dim_, this->nb_dim_);
    sqrtInvHessBarrier(new_sample, new_sqrt_inv_hess_y);

    Dtype scale = r_ / std::sqrt(std::sqrt(Dtype(this->nb_dim_)*Dtype(this->nb_cons_)));
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

  void sqrtInvHessBarrier(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample, Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& new_sqrt_inv_hess){
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> inv_slack = (this->cons_b_ - this->cons_A_ * new_sample).cwiseInverse();

    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> half_hess = inv_slack.asDiagonal()* this->cons_A_;
    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_hess = half_hess.transpose() * half_hess;

    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> new_hess_inv = new_hess.inverse();

    // compute leverage scores
    // Eigen::Matrix<Dtype, Eigen::Dynamic, 1> score  = ((this->cons_A_ * new_hess_inv).cwiseProduct(this->cons_A_)).rowwise().sum().cwiseProduct(inv_slack).cwiseProduct(inv_slack);
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> score  = ((half_hess * new_hess_inv).cwiseProduct(half_hess)).rowwise().sum();

    // compute vaidya hessian
    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> vaidya_new_hess = half_hess.transpose() * score.asDiagonal() * half_hess + Dtype(this->nb_dim_)/Dtype(this->nb_cons_) * new_hess;

    // compute eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> > es(vaidya_new_hess);

    Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic> V = es.eigenvectors();
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> Dv = es.eigenvalues();
    new_sqrt_inv_hess = V * Dv.cwiseInverse().cwiseSqrt().asDiagonal() * V.transpose();
  }

private:
  const Dtype r_;
};

//} // namespace pwalk

#endif // PWALK_VAIDYA_WALKER_HPP_

