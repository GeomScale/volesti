/**
 @file Point.h
 */

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <new>
#include <cmath>
#include <iomanip>
#include <vector>

namespace kdgeraf {
/**
 * A point can be represented by an
 * instance of this class.
 */
template<class DivisionSpace>
class Point {
 public:

  /**
   * The data type we are using.
   */
  typedef typename DivisionSpace::FT FT;

  /**
   * Default constructor.
   */
  Point() {
  }

  /**
   * \brief Constructor to initialize all values with c.
   * Point will be created dynamically,
   * with an array of d size.
   * @param d - dimension of point
   * @param c - value for all the coordinates
   */
  Point(const int d, const FT c)
      : coords(d) {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      coords.at(i) = c;
  }

  /**
   * \brief Constructor to initialize d dimensional point
   * with the values of array c.
   * @param d - dimension of point
   * @param c - array of coordinates
   */
  Point(const int d, const FT c[])
      : coords(d) {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      coords.at(i) = c[i];
  }

  /**
   * \brief Constructor to initialize d dimensional point
   * with the values of vector c.
   * @param d - dimension of point
   * @param c - vector of coordinates
   */
  Point(const int d, const std::vector<FT>& c)
      : coords(d) {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      coords.at(i) = c[i];
  }

  /**
   * \brief Get dimension of point.
   */
  size_t dim() const {
    return coords.size();
  }

  /**
   * \brief Operator that returns the coordinate at the given index. (constant)
   * @param i - index of the coordinate
   * @return the coordinate at index i
   */
  FT operator[](int i) const {
    //assert(i >= 0 && i < (int)coords.size());
    return coords[i];
  }

  /**
   * \brief Operator that returns the coordinate at the given index. (ref)
   * @param i - index of the coordinate
   * @return the coordinate at index i
   */
  FT& operator[](int i) {
    //assert(i >= 0 && i < (int)coords.size());
    return coords[i];
  }

  /**
   * \brief Returns a pointer to the coordinates of this Point.
   * @return the pointer
   */
  const FT* constData() const {
    return &coords[0];
  }

  /**
   * \brief Returns a pointer to the vector of coordinates.
   * @return the reference
   */
  const std::vector<FT>* get_coords() const {
    return &coords;
  }

  /**
   * \brief Minus operator.
   * @param p - the subtrahend
   * @return the difference point
   */
  Point operator -(const Point& p) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] - p[i];
    return Point(D, coordinates);
  }

  /**
   * \brief Plus operator.
   * @param p - the addendum point
   * @return the result point
   */
  Point operator +(const Point& p) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] + p[i];
    return Point(D, coordinates);
  }

  /**
   * \brief Plus operator. (reference)
   * @param a - the addendum
   * @return the result point
   */
  Point operator +(FT& a) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] + a;
    return Point(D, coordinates);
  }

  /**
   * \brief Plus operator.
   * @param a - the addendum
   * @return the result point
   */
  Point operator +(FT a) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] + a;
    return Point(D, coordinates);
  }

  /**
   * \brief Minus operator. (reference)
   * @param a - the subtrahend
   * @return the result point
   */
  Point operator -(FT& a) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] - a;
    return Point(D, coordinates);
  }

  /**
   * \brief Minus operator.
   * @param a - the subtrahend
   * @return the result point
   */
  Point operator -(FT a) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] - a;
    return Point(D, coordinates);
  }

  /**
   * \brief Multiply operator for a single value.
   * All coordinates will be multiplied with the given value.
   * @param m - the multiplier
   * @return the <b> new point </b>
   */
  const Point operator *(const FT m) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] * m;
    return Point(D, coordinates);
  }

  /**
   * \brief Multiply operator for a single value.
   * All coordinates will be multiplied with the given value.
   * @param m - the multiplier
   * @param p - the point to be multiplied
   * @return the new point
   */
  friend const Point operator *(const FT m, const Point& p) {
    return p * m;
  }

  /**
   * \brief Multiply all coordinates of this point with the given value.
   * @param m - the multiplier
   */
  void operator *=(const FT m) {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      coords[i] = coords[i] * m;
  }

  /**
   * \brief Dot product operator.
   * @param p - another point
   * @return the dot product of the two points
   */
  FT operator *(const Point& p) const {
    FT product = 0;
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      product += p[i] * coords[i];
    return product;
  }

  /**
   * \brief Division operator for a single value.
   * All coordinates will be divided by the given value.
   * @param w - the divisor
   * @return the <b> new point </b>
   */
  const Point operator /(const FT w) const {
    const size_t D = coords.size();
    FT coordinates[D];
    for (size_t i = 0; i < D; ++i)
      coordinates[i] = coords[i] / w;
    return Point(D, coordinates);
  }

  /**
   * \brief Division operator for a single value.
   * All coordinates will be divided by the given value.
   * @param w - the divisor
   * @param p - the point to be divided
   * @return the new point
   */
  friend const Point operator /(const FT w, const Point& p) {
    return p / w;
  }

  /**
   * \brief Divide all coordinates of this point by the given value.
   * @param w - the divisor
   */
  void operator /=(const FT w) {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      coords[i] = coords[i] / w;
  }

  /**
   * \brief Equality operator based on the coordinates of the points.
   * @param p - point to compare with
   * @return true if this point is equal to p, else false
   */
  bool operator ==(const Point& p) const {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      if (p[i] != coords[i])
        return false;
    return true;
  }

  /**
   * \brief Returns the norm of this vector.
   * @return the norm
   */
  FT norm() const {
    FT result = 0;
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      result += coords[i] * coords[i];
    //FT r = result;//.to_FT();
    return sqrt(result);
  }

  /**
   * \brief Returns the squared norm of this vector
   * @return the squared norm
   */
  FT squaredNorm() const {
    FT result = 0;
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      result += coords[i] * coords[i];
    return result;
  }

  /**
   * \brief Normalize this point and return a new point with the calculated coordinates.
   * @return the normalized point
   */
  const Point normalized() const {
    return (*this / norm());
  }

  /**
   * \brief Normalize this point.
   */
  void normalize() {
    const FT n = norm();
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      coords[i] = coords[i] / n;
  }

  /** \brief Euclidean distance squared.
   *
   * Faster to compute than Euclidean distance and
   * enough for comparison. This version unrolls
   * the loop.
   *
   * @param p - the other point
   * @return  - the squared Euclidean distance of p
   * and this point
   */
  FT operator &(Point& p) {
    typename std::vector<FT>::iterator it1 = this->coords.begin();
    typename std::vector<FT>::iterator it2 = p.coords.begin();

    FT squared_distance = 0.;
    FT diff0, diff1, diff2, diff3;
    typename std::vector<FT>::iterator last = it1 + this->dim();
    typename std::vector<FT>::iterator lastgroup = last - 3;

    /* Process 4 items with each loop for efficiency. */
    while (it1 < lastgroup) {
      diff0 = it1[0] - it2[0];
      diff1 = it1[1] - it2[1];
      diff2 = it1[2] - it2[2];
      diff3 = it1[3] - it2[3];
      squared_distance += diff0 * diff0 + diff1 * diff1 + diff2 * diff2
          + diff3 * diff3;
      it1 += 4;
      it2 += 4;
    }
    while (it1 < last) {
      diff0 = it1[0] - it2[0];
      squared_distance += diff0 * diff0;
      it1++;
      it2++;
    }

    return squared_distance;
  }

  /**
   * Point iterator.
   */
  typedef typename std::vector<FT>::iterator Point_it;

  /**
   * \brief Iterator pointing at
   * the beginning of coordinates.
   *
   * @return - the iterator
   */
  typename std::vector<FT>::iterator begin() {
    return coords.begin();
  }

  /**
   * \brief Iterator pointing at
   * the end of coordinates.
   *
   * @return - the iterator
   */
  typename std::vector<FT>::iterator end() {
    return coords.end();
  }

  /**
   * \brief Print the point.
   */
  void print() const {
    const size_t D = coords.size();
    for (size_t i = 0; i < D; ++i)
      std::cout << std::setprecision(20) << coords[i] << " ";
  }

  /**
   * \brief Write a point to the given stream.
   * @param s - the stream
   * @param p - the point
   */
  friend std::ostream& operator <<(std::ostream& s, const Point& p) {
    s << std::fixed << std::setprecision(9) << std::setw(12);
    for (size_t i = 0; i < p.dim(); ++i)
      s << p[i] << " ";
    return s;
  }

 private:
  /**
   * The coordinates of the point.
   */
  std::vector<FT> coords;
};
}

#endif  /* POINT_H */
