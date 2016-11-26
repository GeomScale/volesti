/**
 @file Random_generator.h
 */

#ifndef RANDOMGEN_H
#define RANDOMGEN_H

#include <cstdlib>
#include <stdint.h>
#include <random>

/**
 * \brief Generates (uniformly)
 *  a random number in [0, n-1]
 *
 *  @param n - upper bound
 *  @return  - generated number
 */
int myRandomUniform(int n) {
  int top = ((((RAND_MAX - n) + 1) / n) * n - 1) + n;
  int r;
  do {
    r = rand();
  } while (r > top);
  return (r % n);
}


/**
 * \brief Mersenne number generator.
 */
std::mt19937 rng(time(NULL));

/**
 * \brief Generates (uniformly)
 *  a random number in [0, n]
 *
 *  @param n - upper bound
 *  @return  - generated number
 */
int random(int n) {
  std::uniform_int_distribution<int> distribution(0, n);
  return distribution(rng);
}

/**
 * \brief Generates (uniformly)
 *  a random number in [min, max]
 *
 *  @param min - lower bound
 *  @param max - upper bound
 *  @return  - generated number
 */
template<typename T>
int random(T min, T max) {
  std::uniform_real_distribution<T> distribution(min, max);
  return distribution(rng);
}

/**
 * \brief Generates a random point.
 *
 * The random coordinate is generated in [0, n] and
 * has 0.5 possibility to receive a negative sign.
 *
 * @param q - random point
 * @param n - upper bound
 * for random generator
 * @return  - generated number
 */
void random_query(std::vector<float>& q, int n) {
  for(size_t i = 0; i < q.size(); ++i) {
    q[i] = random(n);
    if(myRandomUniform(2))
      q[i] *= -1;
  }
}

#endif /*RANDOMGEN_H*/
