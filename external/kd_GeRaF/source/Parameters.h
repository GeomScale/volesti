/**
 @file Parameters.h
 enum Rotation and Params struct.
 */
 
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

/**
 * Enum for rotation.
 */
enum Rotation {
  /**
   * no rotation
   */
  No,
  /**
   * rotation by using Householder matrices.
   */
  Householder,
};
 
/** @struct Params
 *  @brief Struct of parameters.
 */
typedef struct Params {
  /**
   * max number of points per leaf allowed
   */
  int points_per_leaf;
  /**
   * size of forest
   */
  int trees_no;
  /**
   * number of dimensions to build a tree
   */
  int t;
  /**
   * maximum number of leaves to check
   */
  int max_leaf_check;
  /**
   * flag for rotation
   */
  Rotation rotate_option;
  /**
   * flag for shuffling
   */
  bool shuffle_enable;
  /**
   * for computing variances 
   */
  size_t sample_size;
} 
  /**
   * typedef for struct of Params
   */
Params;

/**
 * Assign the data members of the struct 'parameters'
 * to the corresponding variables.
 *
 * @param parameters      - struct of parameters
 * @param points_per_leaf - max number of points per leaf allowed
 * @param trees_no        - size of forest
 * @param t               - number of dimensions to build a tree
 * @param max_leaf_check  - maximum number of leaves to check
 * @param rotate_option   - flag for rotation
 * @param shuffle_enable  - flag for shuffling
 * @param sample_size     - for computing variances
 */
void assign(const Params& parameters, int& points_per_leaf, int& trees_no, int& t,
            int& max_leaf_check, Rotation& rotate_option, bool& shuffle_enable,
            size_t& sample_size) {
  points_per_leaf = parameters.points_per_leaf;
  trees_no        = parameters.trees_no;
  t               = parameters.t;
  max_leaf_check  = parameters.max_leaf_check;
  rotate_option   = parameters.rotate_option;
  shuffle_enable  = parameters.shuffle_enable;
  sample_size     = parameters.sample_size;
}

#endif /* PARAMETERS_H_ */
