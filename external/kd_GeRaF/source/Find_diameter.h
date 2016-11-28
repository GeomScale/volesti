/**
 @file Find_diameter.h
 This file is NOT used.
 */

#ifndef FIND_DIAMETER_H_
#define FIND_DIAMETER_H_

#include "Division_Euclidean_space.h"


/**
 * \brief Find diameter of a set of points.
 *
 * @param N           - number of points
 * @param D           - dimension of points
 * @param v           - vector of points
 * @return            - diameter
 */
template<typename T>
float find_diameter_exact(const int N, const int D,
                              const std::vector<T>& v) {
  // Note that indexing is of the form: [i * D + j]
  int offset;
  size_t max = 0;
  float d;
  for (int n = 0; n < N; ++n) {
    offset = n * D;
    for (int i = n + 1; i < N; ++i) {
      // d(v[n], v[i])
      d = squared_Eucl_distance<T>(v, offset, offset + D, i * D + D);
      if (d > max) {
        max = d;
      }
    }
  }
  return sqrt(max);
}

/**
 * \brief Find diameter of a set of points.
 *
 * @param v_p         - vector of points
 * @param v_i_begin   - start of points
 * @param v_i_end     - end of points
 * @param D           - dimensions of points
 * @return            - diameter
 */
template<typename T>
float find_diameter_exact(const std::vector<T>& v_p,
                          std::vector<size_t>::iterator& v_i_begin,
                          std::vector<size_t>::iterator& v_i_end,
                          const size_t D) {
  // Note that indexing is of the form: [i * D + j]
  size_t max = 0;
  float d;
  size_t offset1;
  for(std::vector<size_t>::iterator it1 = v_i_begin; it1 != v_i_end; ++it1) {
    offset1 = *it1 * D;
    for(std::vector<size_t>::iterator it2 = it1 + 1; it2 != v_i_end; ++it2) {
      // d(v[i], v[i + 1])
      d = squared_Eucl_distance<T>(v_p, offset1, offset1 + D, *it2 * D);
      if (d > max) {
        max = d;
      }
    }
  }
  return sqrt(max);
}

/**
 * \brief Approximately find diameter of a set of points.
 *
 * Pick random point x. Pick point y farthest from x.
 * Pick point z farthest from y.
 * The distance between y and z can than be used as the
 * diameter.
 *
 * The diameter is no smaller than this value and
 * no larger than twice this value.
 * Proof:
 * Let D = d(p,q) be the diameter. Then d(y,z) ≤ D
 * (since p,q is the maximum argument of d(*,*)),
 * and D = d(p,q) ≤ d(p,y) + d(y,q) ≤ 2d(y,z)
 * (triangle inequality and since z is the maximum
 * argument of d(y,*) = d(*,y)).
 *
 * @param N           - number of points
 * @param D           - dimension of points
 * @param v           - vector of points
 * @return            - diameter
 */
template<typename T>
float find_diameter_appr(const int N, const int D,
                             const std::vector<T>& v) {
  // Note that indexing is of the form: [i * D + j]
  int offset, offset_x = random(N) * D;
  size_t max = 0;
  int point_y_index;
  float d;
  for (int i = 0; i < N; ++i) {
    offset = i * D;
    // d(v[i], x)
    d = squared_Eucl_distance<T>(v, offset, offset + D, offset_x);
    if (d > max) {
      max = d;
      point_y_index = i;
    }
  }
  int offset_y = point_y_index * D;
  max = 0;
  for (int i = 0; i < N; ++i) {
    offset = i * D;
    // d(v[i], y)
    d = squared_Eucl_distance<T>(v, offset, offset + D, offset_y);
    if (d > max) {
      max = d;
    }
  }
  return sqrt(max);
}

/**
 * \brief Approximately find diameter of a set of points.
 *
 * Pick random point x. Pick point y farthest from x.
 * Pick point z farthest from y.
 * The distance between y and z can than be used as the
 * diameter.
 *
 * The diameter is no smaller than this value and
 * no larger than twice this value.
 * Proof:
 * Let D = d(p,q) be the diameter. Then d(y,z) ≤ D
 * (since p,q is the maximum argument of d(*,*)),
 * and D = d(p,q) ≤ d(p,y) + d(y,q) ≤ 2d(y,z)
 * (triangle inequality and since z is the maximum
 * argument of d(y,*) = d(*,y)).
 *
 * @param v_p         - vector of points
 * @param v_i_begin   - start of points
 * @param v_i_end     - end of points
 * @param D           - dimensions of points
 * @return            - diameter
 */
template<typename T>
float find_diameter_appr(const std::vector<T>& v_p,
                              std::vector<size_t>::iterator &v_i_begin,
                              std::vector<size_t>::iterator &v_i_end,
                              const size_t D) {
  // Note that indexing is of vector is of the form: [i * D + j]
  // size = std::distance(itBegin, itEnd) for non-"Random Access" iterators
  size_t N = v_i_end - v_i_begin;
  size_t offset, offset_x = random(N) * D;
  float d, max = 0;
  size_t index_y = 0;
  for(std::vector<size_t>::iterator it = v_i_begin; it != v_i_end; ++it) {
    offset = *it * D;
    d = squared_Eucl_distance<T>(v_p, offset, offset + D, offset_x);
    if(d > max) {
      max = d;
      index_y = *it;
    }
  }
  max = 0;
  size_t offset_y = index_y * D;
  for(std::vector<size_t>::iterator it = v_i_begin; it != v_i_end; ++it) {
    offset = *it * D;
    d = squared_Eucl_distance<T>(v_p, offset, offset + D, offset_y);
    if(d > max) {
      max = d;
    }
  }
  return sqrt(max);
}

#endif /* FIND_DIAMETER_H_ */
