/**
 @file Mean_variance.h
 */


/*
 * File:   Mean_variance.h
 * Author: SAMARAS
 *
 * Created on 2014
 */

#ifndef MEAN_VARIANCE_H
#define MEAN_VARIANCE_H

#include <cmath>
#include "Find_k_max.h"
#include "Division_Euclidean_space.h"
#include <cstring>
#include <cassert>

/**
 * \brief Calculate the average and scaled variance along each dimension
 * and return the 't' coordinates with highest variances. Keeps the average
 * and variance arrays at local scope.
 *
 * @param points  - vector of points
 * @param size    - number of points
 * @param D       - dimension of points
 * @param var  		- array of variances (to be filled by this function)
 * @param sample  - take into account only this sample from the whole dataset,
 * value 100 seems to work well. If -1, then the whole dataset is taken into account.
 */
template<typename T>
void compute_variances(const std::vector<T>& points, size_t size, const int D,
                       float* var, int sample = -1)
{
  assert(points.size() > 0);
  if(sample < 0)
    sample = size;
  else if((int)size > sample)
    size = sample;
    
	float avg[D];

  // i = 0, i_n = 1
  assert(D > 0);
  #if __cplusplus >= 201103L
    std::copy_n(&points[0], D, avg);
  #else
    std::copy(&points[0], &points[0] + D, avg);
  #endif

  // i = 1, i_n = 0.5
  if (size >= 2) {
    //assert(points[1].dim() == D);
    for (int d = D - 1; d >= 0; --d) {
      float const delta = points[1 * D + d] - avg[d];
      avg[d] += delta * 0.5f;
      var[d] = delta * (points[1 * D + d] - avg[d]);
    }
  } else {
    std::fill_n(var, D, 0.0f);
  }

  // i = 2, ...
  for (size_t i = 2; i < size; ) {
    {
      const float i_n = 1.0f / (1.0f + i);
      //assert(points[i].dim() == D);
      for (int d = 0; d < D; ++d) {
        float const delta = points[i * D + d] - avg[d];
        avg[d] += delta * i_n;
        var[d] += delta * (points[i * D + d] - avg[d]);
      }
    }
    ++i;

    if (i >= size) break;
    {
      const float i_n = 1.0f / (1.0f + i);
      //assert(points[i].dim() == D);
      for (int d = D - 1; d >= 0; --d) {
        float const delta = points[i * D + d] - avg[d];
        avg[d] += delta * i_n;
        var[d] += delta * (points[i * D + d] - avg[d]);
      }
    }
    ++i;
  }
}

/**
 * \brief Calculate the average and scaled variance along each dimension
 * and return the 't' coordinates with highest variances.
 * @param t           - number of highest variances
 * @param points      - vector of points
 * @param D           - dimension of points
 * @param size        - number of points
 * @param avg         - array of averages
 * @param var         - array of variances
 * @param split_dims  - the t dimensions with the highest variances (output)
 * @param sample      - take into account only this sample from the whole dataset,
 * value 100 seems to work well. If -1, then the whole dataset is taken into account.
 */
template<typename T>
void compute_variances(size_t t, const std::vector<T>& points, const int D,
                       size_t size, float* avg,
                       float* var, size_t* split_dims, int sample = -1)
{
  assert(points.size() > 0);
  if(sample < 0)
    sample = size;
  else if((int)size > sample)
    size = sample;

  // i = 0, i_n = 1
  assert(D > 0);
  #if __cplusplus >= 201103L
    std::copy_n(&points[0], D, avg);
  #else
    std::copy(&points[0], &points[0] + D, avg);
  #endif

  // i = 1, i_n = 0.5
  if (size >= 2) {
    //assert(points[1].dim() == D);
    for (int d = D - 1; d >= 0; --d) {
      float const delta = points[1 * D + d] - avg[d];
      avg[d] += delta * 0.5f;
      var[d] = delta * (points[1 * D + d] - avg[d]);
    }
  } else {
    std::fill_n(var, D, 0.0f);
  }

  // i = 2, ...
  for (size_t i = 2; i < size; ) {
    {
      const float i_n = 1.0f / (1.0f + i);
      //assert(points[i].dim() == D);
      for (int d = 0; d < D; ++d) {
        float const delta = points[i * D + d] - avg[d];
        avg[d] += delta * i_n;
        var[d] += delta * (points[i * D + d] - avg[d]);
      }
    }
    ++i;

    if (i >= size) break;
    {
      const float i_n = 1.0f / (1.0f + i);
      //assert(points[i].dim() == D);
      for (int d = D - 1; d >= 0; --d) {
        float const delta = points[i * D + d] - avg[d];
        avg[d] += delta * i_n;
        var[d] += delta * (points[i * D + d] - avg[d]);
      }
    }
    ++i;
  }

  /* Find t dimensions with largest scaled variance. */
  kthLargest(var, D, t, split_dims);
  //for(int g=0;g<t;++g)
  //	std::cout << var[g] << "\n";
 
  /*
  klein
  2703.59
2685.96
3596.94
3575.21
3578.98
sphere
557.436
564.388
566.614
563.327
575.483

*/
  	
}

/**
 * \brief Find the median element of a collection of points,
 * in a specific dimension.
 *
 * Example:
 * Input: the indices:
 * 0 1 2 3 4 5 6 7
 *
 * Output:
 * median = 4
 * 2 1 0 3 4 5 6 7
 *
 * with this dataset:
 * 0 0
 * 1 1
 * 2 2
 * 3 3
 * 4 4
 * 5 5
 * 6 6
 * 7 7
 *
 * @param v_p       - data set
 * @param v_i_begin - begin of vector of indices
 * @param v_i_end   - end of vector of indices
 * @param middle    - index of the middle element
 * @param atDim     - dimension in which to search
 * @param D         - dimension of data set
 * @return          - the median
 */
template<typename T>
std::vector<size_t>::iterator find_median(
    const std::vector<T>& v_p, std::vector<size_t>::iterator &v_i_begin,
    std::vector<size_t>::iterator &v_i_end, const size_t middle,
    const size_t atDim, const size_t D)
{
  std::nth_element(v_i_begin, v_i_begin + middle, v_i_end,
                   [&](int lhs, int rhs) {
                     return v_p[lhs * D + atDim] < v_p[rhs * D + atDim];
                   });
  //std::cout << "middle = " << middle << std::endl;
  return v_i_begin + middle;
}

/**
 * \brief Find the median element of a collection of points,
 * in a specific dimension.
 *
 * @param v_i_begin - begin of vector of indices
 * @param v_i_end   - end of vector of indices
 * @param middle    - index of the middle element
 * @param t_dims    - specifies the dimension
 * @return          - the median
 */
template<typename T>
std::vector<size_t>::iterator find_median(
    std::vector<size_t>::iterator &v_i_begin,
    std::vector<size_t>::iterator &v_i_end, const size_t middle,
    std::vector<T>& t_dims)
{
  std::nth_element(v_i_begin, v_i_begin + middle, v_i_end,
                   [&t_dims](int lhs, int rhs) {
                     return t_dims[lhs] < t_dims[rhs];
                   });
  return v_i_begin + middle;
}

#endif  /* MEAN_VARIANCE_H */
