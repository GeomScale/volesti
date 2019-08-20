#ifndef PWALK_WALK_HPP_
#define PWALK_WALK_HPP_

#include <Eigen/Dense>
//#include "common.h"

//namespace pwalk {

template <typename Dtype>
class Walker {
public:
  Walker(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& initialization, const Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& cons_A, const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& cons_b) : nb_dim_(cons_A.cols()), nb_cons_(cons_A.rows()), nb_curr_samples_(1), initialization_(initialization), cons_A_(cons_A), cons_b_(cons_b), curr_sample_(initialization){}

  virtual ~Walker(){}

  virtual bool doSample(Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample, Dtype lazy = Dtype(0.5)){
    return false;
  }

  // check whether a given point is in th polytope
  bool checkInPolytope(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& new_sample){
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
  Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& getCurrSample() {
    return curr_sample_;
  }

protected:

  // Dimension
  const int nb_dim_;
  // number of constraints
  const int nb_cons_;
  // current sample size
  int nb_curr_samples_;
  // Initial vector
  const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& initialization_;
  // constraints matrix A
  const Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& cons_A_;
  // constraint vector b
  const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& cons_b_;
  // Current vector
  Eigen::Matrix<Dtype, Eigen::Dynamic, 1> curr_sample_;
};


//} // namespace pwalk

#endif // PWALK_WALK_HPP_
