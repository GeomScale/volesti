/**
 @file Division_Euclidean_space.h
 A space that allows us to split the data.
 */

#ifndef DIVISION_EUCLIDEAN_SPACE_H
#define DIVISION_EUCLIDEAN_SPACE_H

#include <iostream>

/**
 * \implements DivisionSpace
 *
 * Division Euclidean space hosts
 * the points of the data set.
 *
 * Internally, the points are kept
 * in a 1D vector and in order to
 * access the j-th coordinate
 * of the i-th point, we do:
 *
 * p[i * D + j]
 */
template<typename T>
class Division_Euclidean_space {
 public:
  /**
   * The data type.
   */
  typedef T FT;

  /**
   * Constructor, which
   * sets 'N' and 'D' to zero.
   */
  Division_Euclidean_space()
      : N(0),
        D(0) {

  }

  /**
   * Constructor, which
   * reserves space for the data.
   *
   * @param n - size of the data
   * @param d - dimension of the data
   */
  Division_Euclidean_space(const size_t& n, int d)
      : N(n),
        D(d) {
    p.reserve(N * D);
  }

  /**
   * @param n - size of data
   */
  void setSize(size_t& n) {
    N = n;
  }
  
  /**
   * @param n - size of data
   */
  void setSize(int n) {
    N = n;
  }

  /**
   * @param d - dimension of data
   */
  void setDim(int& d) {
    D = d;
  }

  /**
   * \brief Inserts a new value to the collection of
   * points, held in the private vector.
   *
   * @param v - value to be inserted
   */
  void insert(FT v) {
    try {
  		p.push_back(v);
		} catch (const std::bad_alloc&) {
			std::cerr << "Out of memory when reading data, exiting..." << std::endl;
  		exit(1);
		}
  }

  /**
   * Get the vector of points
   *
   * @return - the vector
   */
  const std::vector<FT>& getPoints() const {
    return p;
  }
  
  /**
   * Get the vector of the sum of squared
   * coordinates per point.
   *
   * @return - the vector
   */
  const std::vector<FT>& getSquaredCoords() const {
    return squared_coords;
  }

  /**
   * Get the vector of points
   *
   * @return - the vector
   */
  std::vector<FT>* getPointsPointer() {
    return &p;
  }

  /**
   * Get the number of points
   *
   * @return - the number of points
   */
  const size_t& size() const {
    return N;
  }

  /**
   * Get the dimension of points (ref.)
   *
   * @return - the dimension
   */
  const int& dim() const {
    return D;
  }

  /**
   * Get the dimension of points
   *
   * @return - the dimension
   */
  const int dim() {
    return D;
  }

  /**
   * Print the points.
   */
  void printPoints() {
    int c = 0;
    for (size_t i = 0; i < p.size(); ++i) {
      std::cout << p[i] << " ";
      if (c++ == D - 1) {
        std::cout << "\n";
        c = 0;
      }
    }
  }

  /**
   * Compute and store the sum of squared coordinates per point.
   */  
  void computeSquare() {
    int c = 0;
    FT sq = 0;
    squared_coords.resize(N);
    for (size_t i = 0; i < p.size(); ++i) {
      sq += p[i];
      if (c++ == D - 1) {
        squared_coords.push_back(sq);
        sq = 0;
        c = 0;
      }
    }
  }  

 private:
  /**
   * number of points
   */
  size_t N;
  /**
   * dimension of points
   */
  int D;
  /**
   * vector of points
   * Note that indexing is of the form: [i * D + j]
   */
  std::vector<FT> p;
  /**
   * vector of sum of the squared coordinates of every point
   * At index 0 lies the sum of the squared coordinates of the
   * first point and so on.
   */
  std::vector<FT> squared_coords;
};

////////////////////////////////////////////////////////
////////////////// END OF CLASS ////////////////////////
////////////////////////////////////////////////////////

/** \brief Euclidean distance squared.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison. Parameters are vectors.
 *
 * @param p1 - first point
 * @param p2 - second point
 * @return - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distance(std::vector<T> &p1,
                                std::vector<T> &p2) {
  typename std::vector<T>::iterator it2 = p2.begin();
  float squared_distance = 0.;
  float diff;
  for (typename std::vector<T>::iterator it1 = p1.begin(); it1 < p1.end();
      ++it1, ++it2) {
    diff = *it1 - *it2;
    squared_distance += diff * diff;
  }
  return squared_distance;
}

/* TODO: Make all searches in the tree use squared Eucl distance. */

/** \brief Euclidean distance squared.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison. Parameters are vectors.
 * First vector is a collection of points, thus
 * we need its start and end.
 *
 * @param p1    - first point
 * @param start - starting index of first point
 * @param end  - ending index of first point
 * @param p2    - second point
 * @return      - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distance(const std::vector<T> &p1, size_t start,
                                size_t end, const std::vector<T> &p2) {
               // this one
  typename std::vector<T>::const_iterator it2 = p2.begin();
  T squared_distance = 0;
  T diff;
  for (size_t i = start; i < end; ++i, ++it2) {
    diff = p1[i] - *it2;
    squared_distance += diff * diff;
  }
  return squared_distance;
}

/** \brief Euclidean distance squared, unroll the loop.
 * This is specialized, avoid usage.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison. Parameters are vectors.
 * First vector is a collection of points, thus
 * we need its start and end.
 *
 * @param p1    - first point
 * @param start - starting index of first point
 * @param end  - ending index of first point
 * @param p2    - second point
 * @return      - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distanceU(const std::vector<T> &p1, size_t start,
                                size_t end, const std::vector<T> &p2) {
/* meant to be used only for 128 dimensions, but it is still slower.
First  method: 2.49e-07 seconds.
Second method: 7.88e-07 seconds.

First  method: 2.51e-07 seconds.
Second method: 7.89e-07 seconds.

First  method: 2.55e-07 seconds.
Second method: 6.65e-07 seconds.
*/

  typename std::vector<T>::const_iterator it2 = p2.begin();
  T squared_distance = 0;
  T diff0, diff1, diff2, diff3;
  for (size_t i = start; i < end; ++i, ++it2) {
    diff0 = p1[i++] - *it2++;
    diff1 = p1[i++] - *it2++;
    diff2 = p1[i++] - *it2++;
    diff3 = p1[i] - *it2;
    
    squared_distance += (diff0 * diff0) + (diff1 * diff1) + (diff2 * diff2) + (diff3 * diff3);
  }
  return squared_distance;
}

/** \brief Euclidean distance squared.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison. Parameters are vectors.
 * We need to compute only the dot product.
 *
 * @param p1    - first point
 * @param start - starting index of first point
 * @param end  - ending index of first point
 * @param p1_sq - squared first point
 * @param p2    - second point
 * @param p2_sq - squared second point
 * @return      - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distance(const std::vector<T> &p1, size_t start, size_t end,
													T p1_sq, const std::vector<T> &p2, T p2_sq) {
/*
Slower than the usual method:
First  method: 1.52e-07 seconds.
Second method: 4.89e-07 seconds.

First  method: 1.54e-07 seconds.
Second method: 5.03e-07 seconds.
*/
  typename std::vector<T>::const_iterator it2 = p2.begin();
  T dot_product = 0;
  for (size_t i = start; i < end; ++i, ++it2) {
  	dot_product += p1[i] * *it2;
  }
  return p1_sq + p2_sq - 2 * dot_product;
}

/** \brief Euclidean distance squared.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison.
 *
 * We pass a vector of points, where both points
 * lie into, thus we need the start and the end
 * of of the first point and the start of the
 * second one, since points should be of same
 * dimension.
 *
 * @param p       - vector of points
 * @param start1  - starting index of first point
 * @param end1   - ending index of first point
 * @param start2  - starting index of second point
 * @return        - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distance(const std::vector<T> &p, size_t start1,
                                size_t end1, size_t start2) {
  typename std::vector<T>::const_iterator it2 = p.begin() + start2;
  float squared_distance = 0.;
  float diff;
  for (size_t i = start1; i < end1; ++i, ++it2) {
    diff = p[i] - *it2;
    squared_distance += diff * diff;
  }
  return squared_distance;
}

/** \brief Euclidean distance squared.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison. Parameters are iterators.
 *
 * We pass a vector of points, where both points
 * lie into, thus we need the start and the end
 * of of the first point and the start of the
 * second one, since points should be of same
 * dimension.
 *
 * @param p       - vector of points
 * @param start1  - starting of first point
 * @param end1   - ending of first point
 * @param start2  - starting of second point
 * @return        - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distance(const std::vector<T> &p,
                        std::vector<size_t>::iterator& start1,
                        std::vector<size_t>::iterator end1,
                        std::vector<size_t>::iterator start2) {
  float squared_distance = 0.;
  float diff;
  for (std::vector<size_t>::iterator it = start1; it < end1; ++it, ++start2) {
    diff = *it - *start2;
    squared_distance += diff * diff;
  }
  return squared_distance;
}

/** \brief Euclidean distance squared.
 *
 * Faster to compute than Euclidean distance and
 * enough for comparison. Parameters are iterators
 * for the first point and index for the second one.
 *
 * We pass a vector of points, where both points
 * lie into, thus we need the start and the end
 * of of the first point and the index of the
 * second one, since points should be of same
 * dimension.
 *
 * @param p       - vector of points
 * @param start1  - starting of first point
 * @param end1   - ending of first point
 * @param index2  - index of second point
 * @return        - the Euclidean distance of p1-p2
 */
template<typename T>
T squared_Eucl_distance(const std::vector<T> &p,
                        std::vector<size_t>::iterator& start1,
                        std::vector<size_t>::iterator end1,
                        size_t index2) {
  float squared_distance = 0.;
  float diff;
  for (std::vector<size_t>::iterator it = start1; it < end1; ++it, ++index2) {
    diff = *it - p[index2];
    squared_distance += diff * diff;
  }
  return squared_distance;
}

/** \brief Euclidean distance.
 *  @param p1 - first point
 *  @param start - starting index of first point
 *  @param end  - ending index of first point
 *  @param p2 - second point
 *  @return - the Euclidean distance
 */
template<typename T>
T Eucl_distance(const std::vector<T> &p1, size_t start,
                        size_t end, const std::vector<T> &p2) {
  return std::sqrt(squared_Eucl_distance(p1, start, end, p2));
}

/**
 * \brief Compute and store the sum of squared coordinates per point.
 * @param p								-	input vector of points
 * @param N 							- number of points
 * @param D 							- dimension
 * @param squared_coords 	- output vector
 */
template<typename T>
void computeSquare(const std::vector<T>& p, const size_t N, const size_t D,
                   std::vector<T>& squared_coords) {
  unsigned int c = 0;
  T sq = 0;
  for (size_t i = 0; i < p.size(); ++i) {
    sq += p[i] * p[i];
    if (c++ == D - 1) {
      squared_coords.push_back(sq);
      sq = 0;
      c = 0;
    }
  }
}

/**
 * \brief Compute and store the sum of squared coordinates per point.
 * Here, the input vector is 2D.
 *
 * @param p								-	input vector of points
 * @param N 							- number of points
 * @param D 							- dimension
 * @param squared_coords 	- output vector
 */
template<typename T>
void computeSquare(const std::vector< std::vector<T> >& p, const size_t N, const size_t D,
                   std::vector<T>& squared_coords) {
  T sq = 0;
  for (size_t i = 0; i < p.size(); ++i) {
    for (size_t j = 0; j < D; ++j) {
      sq += p[i][j] * p[i][j];
    }
    squared_coords.push_back(sq);
    sq = 0;
  }
}

#endif /*DIVISION_EUCLIDEAN_SPACE_H*/
