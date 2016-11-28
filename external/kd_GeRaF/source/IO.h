/**
 @file IO.h
 */

#ifndef IO_H
#define IO_H

#include <fstream>
#include "Point.h"
using namespace kdgeraf;

/** \brief Read an Euclidean division space from file.
 *
 * First line of file should have 'dimension number_of_points'.
 *
 * @param ds       - Euclidean division space
 * @param filename - input file
 */
template<typename T>
void readDivisionSpace(Division_Euclidean_space<T>& ds, char* filename) {
  std::ifstream infile;
  int N = 0;
  int D;

  infile.open(filename);
  if (!infile)
    std::cout << "File not found!" << std::endl;
  // Read first line of file
  infile >> D;
  if (D < 1) {
    std::cout << "ERROR, dimension less than one!!\n\n";
    return;
  }
  infile >> N;
  if (N < 1) {
    std::cout << "ERROR, number of points less than 1!!\n\n";
    return;
  }

  ds.setSize((size_t&)N);
  ds.setDim(D);

  int hRead = 0;
  float coord;
  for (int n = 0; n < N && infile; ++n) {
    for (int i = 0; i < D; ++i) {
      infile >> coord;
      ds.insert(coord);
    }
    hRead++;
  }
  if (hRead != N)
    std::cout << "ERROR, read less than " << N << " points!!\n\n";
}

/** \brief Read an Euclidean division space from file.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param ds       - Euclidean division space
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void readDivisionSpace(Division_Euclidean_space<T>& ds, int N, int D, char* filename) {
  std::ifstream infile;

  infile.open(filename);
  if (!infile)
    std::cout << "File not found!" << std::endl;

  int hRead = 0;
  float coord;
  for (int n = 0; n < N && infile; ++n) {
    for (int i = 0; i < D; ++i) {
      infile >> coord;
      ds.insert(coord);
    }
    hRead++;
  }
  if (hRead != N)
    std::cout << "ERROR, read less than " << N << " points!!\n\n";
}

/** \brief Read an Euclidean division space from a "fvecs" file.
 * Every line of the file has this format: D x_1 ... x_D
 *
 * Dimension and number of points are computed by the function.
 *
 * @param ds       - Euclidean division space
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void readDivisionSpacefvecs(Division_Euclidean_space<T>& ds, int& N, int& D, char* filename) {
  FILE* fid;
  fid = fopen(filename, "rb");

  if (!fid) {
    printf("I/O error : Unable to open the file %s\n", filename);
    std::cerr << "Error: " << strerror(errno) << std::endl;
  }

  // we assign the return value of fread() to 'sz' just to suppress a warning
  size_t sz = fread(&D, sizeof(D), 1, fid);
  fseek(fid, 0L, SEEK_END);
  sz = ftell(fid);
  N = sz / (1 * 4 + D * 4);
  //printf("N = %d, D = %d, |%s|\n", N, D, filename);

  fseek(fid, 0L, SEEK_SET);
  ds.setSize(N);
  ds.setDim(D);
  //std::cout << ds.dim() << " " << ds.size() << "\n";
  int c = 0;
  float v;
  int i, j;
  for(i = 0; i < N; ++i) {
    sz = fread(&D, sizeof(D), 1, fid);
    //printf("%d\n", D);
    for (j = 0; j < D; ++j) {
      sz = fread(&v, sizeof(v), 1, fid);
      //if(c >= 279619)
      //printf("j = %d, v = %f, read up to point %d\n", j, v, c);
      	ds.insert(v,c);
    }
    ++c;
    //printf("read up to %d\n", c);
  }
  if(c != N)
  	printf("WARNING! Read less points than expected.\n");
}

/*
j = 255, v = 0.052300, read up to point 279620
j = 256, v = 0.052300, read up to point 279620
terminate called after throwing an instance of 'std::bad_alloc'
  what(): std::bad_alloc
Aborted
Dataset: 500,000 points in 960 dimensions.
*/

/** \brief Read a collection of points from file.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 * Points lie in an Euclidean division space.
 *
 * @param v        - vector of points
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void read_points(std::vector< Point<Division_Euclidean_space<T>> >& v, int N, int D, const char* filename) {
  std::ifstream infile;

  infile.open(filename);
  if (!infile)
    std::cout << "File not found!" << std::endl;

  int hRead = 0;
  T coords[D];
  for (int n = 0; n < N && infile; ++n) {
    for (int i = 0; i < D; ++i)
      infile >> coords[i];
    Point<Division_Euclidean_space<T> > p(D, coords);
    v.push_back(p);
    hRead++;
  }
  if (hRead != N)
    std::cout << "ERROR, read less than " << N << " points!!\n\n";
}

/** \brief Read a collection of points from file.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param v        - vector of points
 * @param N        - number of points
 * @param D        - dimension of points
 * @param filename - input file
 */
template<typename T>
void read_points(std::vector< std::vector<T> >& v, int N, int D, const char* filename) {
  std::ifstream infile;

  infile.open(filename);
  if (!infile)
    std::cout << "File not found!" << std::endl;

  v.resize(N);
  int hRead = 0;
  T coord;
  for (int n = 0; n < N && infile; ++n) {
    for (int i = 0; i < D; ++i) {
      infile >> coord;
      v[n].push_back(coord);
    }
    hRead++;
  }
  if (hRead != N)
    std::cout << "ERROR, read less than " << N << " points!!\n\n";
}

/** \brief Read a collection of indices of points from file.
 *
 * Dimension and number of points should have been
 * assigned a value before reaching this function.
 *
 * @param v        - vector of indices
 * @param N        - number of indices
 * @param k        - dimension of indices
 * @param filename - input file
 */
void read_indices(std::vector<std::vector<int> >& v, int N, int k, const char* filename) {
  std::ifstream infile;

  infile.open(filename);
  if (!infile)
    std::cout << "File not found!" << std::endl;

  v.resize(N);
  int hRead = 0;
  int index;
  for (int n = 0; n < N && infile; ++n) {
    for (int i = 0; i < k; ++i) {
      infile >> index;
      v[n].push_back(index);
    }
    hRead++;
  }
  if (hRead != N)
    std::cout << "ERROR, read less than " << N << " indices!!\n\n";
}

/** \brief Print 1D vector.
 *
 * @param v        - vector to be printed
 */
template<typename T>
void print_vector(std::vector<T>& v) {
  for(int i = 0; i < v.size(); ++i) {
    std::cout << v[i] << std::endl;
  }
}

#endif /*IO_H*/
