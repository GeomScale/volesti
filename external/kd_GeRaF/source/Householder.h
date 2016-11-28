/**
 @file Householder.h
 Performs random rotation.
 */

#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "Division_Euclidean_space.h"
#include "Random_generator.h"

/**
 * Xorshift128 PRNG
 */
typedef struct {
  /**
   * Data of the state.
   */
  uint64_t u64[2];
} prng_state;

/**
 * \brief Initialize PRNG.
 *
 * @param p - PRNG state
 */
void prng_init(prng_state * const p) {
  p->u64[0] = (uint64_t) 0x082cdc0b3649c1b8ULL;
  p->u64[1] = (uint64_t) 0x19a7a27c89201352ULL;
}

/**
 * \brief Perturb the PRNG state.
 *
 * @param p - PRNG state
 * @param v - perturbation parameter
 */
static inline void prng_perturb(prng_state * const p, const uint64_t v) {
  p->u64[0] ^= v ^ (v << 11U) ^ (v >> 7U) ^ (v << 23U);
  p->u64[1] ^= (v >> 17U) ^ (v << 31U);
}

/**
 * Next PRNG state.
 *
 * @param p - PRNG state
 * @return  - number helpful for 'prng_delta'
 */
static inline uint64_t prng_next(prng_state * const p) {
  const uint64_t tmp0 = p->u64[1];
  const uint64_t tmp1 = p->u64[0] ^ (p->u64[0] << 23U);
  p->u64[0] = tmp0;
  return (p->u64[1] = (tmp1 ^ tmp0 ^ (tmp1 >> 17U) ^ (tmp0 >> 26U))) + tmp0;
}

/**
 * @param p - PRNG state
 * @return  - PRNG delta
 */
static inline float prng_delta(prng_state * const p) {
  float f;
  do {
    f = (float) prng_next(p) / 9223372036854775808.0f - 1.0f;
  } while (f <= -1.0f || f >= +1.0f);
  return f;
}

/** \brief Generate a pseudorandom unit vector
 * @param p - PRNG state to use
 * @param v - std::vector to save the vector to
 * @param D - Number of dimensions in the vector
 * @return - 0 if success, -1 if PRNG state is stuck on zero.
 */
int householder_vector(prng_state * const p, std::vector<float>& v,
                       const int D) {
  while (1) {
    float norm = 0.0f;
    float scale;

    v.resize(D);

    /* Set components, and calculate 2-norm. */
    for (int d = 0; d < D; ++d) {
      const float c = prng_delta(p);
      v[d] = c;
      norm += c * c;
    }

    /* Xorshift128 failure mode is all zeros. */
    if (norm <= 0.0f)
      return -1;

    /* Calculate the scaling factor, */
    scale = sqrt(1.0 / norm);
    /* scale the components, and recalculate 2-norm. */
    norm = 0.0f;
    for (int d = 0; d < D; ++d) {
      const float c = scale * v[d];
      v[d] = c;
      norm += c * c;
    }

    /* If the 2-norm is correct to 16 bits or so, accept it. */
    if (norm > 0.999984741211f && norm < 1.00001525879f)
      return 0;
  }
  return -1;
}

/** \brief Transform one vector
 * @param N - number of points
 * @param D - Number of dimensions
 * @param transformed - Transformed vector
 * @param original -  Original vector
 * @param vector -  Householder transformation
 * vector (reflection hyperplane unit normal).
 */
template<typename T>
void householder_transform(const int N, const int D,
                           std::vector<T>& transformed,
                           const std::vector<T>& original,
                           const std::vector<float>& vector) {
  for (int i = 0; i < N; ++i) {
    const int offset = i * D;
    float distance = 0.0f;

    /* Dot product between Householder vector and original point
     * == distance from point to Householder reflection hyperplane. */
    for (int d = 0; d < D; ++d)
      distance += vector[d] * original[offset + d];

    /* Transformed point is on the other side of the hyperplane.
     * Double the distance. */
    distance += distance;

    /* Calculate transformed coordinates. */
    for (int d = 0; d < D; ++d)
      transformed.push_back(original[offset + d] - distance * vector[d]);
  }
}

/** \brief Transform one dataset
 * @param N - number of points
 * @param D - Number of dimensions
 * @param transformed - Transformed dataset
 * @param original -  Original dataset
 * @param vector -  Householder transformation
 * vector (reflection hyperplane unit normal).
 */
template<typename T>
void householder_transform(const int N, const int D,
                           Division_Euclidean_space<T>& transformed,
                           const std::vector<T>& original,
                           const std::vector<float>& vector) {
  for (int i = 0; i < N; ++i) {
    const int offset = i * D;
    float distance = 0.0f;

    /* Dot product between Householder vector and original point
     * == distance from point to Householder reflection hyperplane. */
    for (int d = 0; d < D; ++d)
      distance += vector[d] * original[offset + d];

    /* Transformed point is on the other side of the hyperplane.
     * Double the distance. */
    distance += distance;

    /* Calculate transformed coordinates. */
    for (int d = 0; d < D; ++d)
      transformed.insert(original[offset + d] - distance * vector[d]);
  }
}

/** \brief Transform one dataset
 * @param N - number of points
 * @param D - Number of dimensions
 * @param transformed - Transformed dataset
 * @param original -  Original dataset
 * @param vector -  Householder transformation
 * vector (reflection hyperplane unit normal).
 * @param split_dims - transform the dataset only in these dimensions
 * @param t - size of split_dims parameter
 */
template<typename T>
void householder_transform(const int N, const int D,
                           Division_Euclidean_space<T>& transformed,
                           const std::vector<T>& original,
                           const std::vector<float>& vector,
                           size_t* split_dims, const int t) {
  for (int i = 0; i < N; ++i) {
    const int offset = i * D;
    //const int offset_t = i * t;
    float distance = 0.0f;

    /* Dot product between Householder vector and original point
     * == distance from point to Householder reflection hyperplane. */
    for (int d = 0; d < D; ++d)
      distance += vector[d] * original[offset + d];

    /* Transformed point is on the other side of the hyperplane.
     * Double the distance. */
    distance += distance;

    /* Calculate transformed coordinates. */
    for (int d = 0; d < t; ++d)
      //transformed[offset_t + d]
      transformed.insert(
          original[offset + split_dims[d]] - distance * vector[split_dims[d]]);
  }
}

/**
 * \brief Rotate a set of points.
 *
 * @param N           - number of points
 * @param D           - dimension of points
 * @param transformed - transformed/rotated points
 * @param original    - original points
 * @param v           - reflection vector
 * @return            - error flag
 */
template<typename T>
int rotate(const int N, const int D, std::vector<T>& transformed,
           const std::vector<T>& original, std::vector<float>& v) {
  // this is used in auto-configuration
  prng_state prng;
  prng_init(&prng);
  prng_perturb(&prng, myRandomUniform(4000));
  prng_perturb(&prng, time(NULL));

  if (householder_vector(&prng, v, D)) {
    std::cout << "Xorshift128 failure mode is all zeros.\n";
    return -1;
  }
  householder_transform<T>(N, D, transformed, original, v);
  return 0;
}

/**
 * \brief Rotate a set of points.
 *
 * @param N           - number of points
 * @param D           - dimension of points
 * @param transformed - Space where
 * the resulted points will lie
 * @param original    - original points
 * @param v           - reflection vector
 * @return            - error flag
 */
template<typename T>
int rotate(const int N, const int D, Division_Euclidean_space<T>& transformed,
           const std::vector<T>& original, std::vector<float>& v) {
  prng_state prng;
  prng_init(&prng);
  prng_perturb(&prng, myRandomUniform(4000));
  prng_perturb(&prng, time(NULL));

  if (householder_vector(&prng, v, D)) {
    std::cout << "Xorshift128 failure mode is all zeros.\n";
    return -1;
  }
  householder_transform<T>(N, D, transformed, original, v);
  return 0;
}

/**
 * \brief Rotate a set of points.
 *
 * @param N           - number of points
 * @param D           - dimension of points
 * @param transformed - Space where
 * the resulted points will lie
 * @param original    - original points
 * @param split_dims  - transform the data set only in these dimensions
 * @param t           - size of split_dims parameter
 * @param v           - reflection vector
 * @return            - error flag
 */
template<typename T>
int rotate(const int N, const int D, Division_Euclidean_space<T>& transformed,
           const std::vector<T>& original, size_t* split_dims,
           const int t, std::vector<float>& v) {
  prng_state prng;
  prng_init(&prng);
  prng_perturb(&prng, myRandomUniform(4000));
  prng_perturb(&prng, time(NULL));

  if (householder_vector(&prng, v, D)) {
    std::cout << "Xorshift128 failure mode is all zeros.\n";
    return -1;
  }

  householder_transform<T>(N, D, transformed, original, v, split_dims, t);
  return 0;
}

#endif /* HOUSEHOLDER_H */
