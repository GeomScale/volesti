/**
 @file Auto_random_kd_forest.h
 */

#ifndef AUTO_RANDOM_KD_TREE_H
#define AUTO_RANDOM_KD_TREE_H

#include <algorithm>
#include <tuple>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <ctime>


#ifdef RKD_PAR
#include <thread>         // std::thread
#if !(__GNUC__ >= 5 || __GNUC_MINOR__ >= 8)
#include "Hardware_con.h"
#endif
#endif

#include "Householder.h"
#include "Tree.h"
#include "IO.h"
#include "Point.h"
#include "Parameters.h"

using namespace kdgeraf;
/**
 * \brief Compare function for sorting, based
 * on first key.
 *
 * @param a - first pair
 * @param b - second pair
 * @return  - true if first pair
 * is less than the second. Otherwise, false.
 */
bool lib_compar(const std::pair<float, int>& a,
                const std::pair<float, int>& b) {
  return (a.first < b.first);

}

/**
 * A forest of Tree's. The points of the trees,
 * lie in the DivisionSpace. This class receives
 * the dataset and the queries, builds the forest
 * automatically (by performing automatic computation
 * of the optimal parameters), performs the search and
 * finally returns the results.
 */
template<typename T>
class Auto_random_kd_forest {
 private:

  /** \brief Read a data set, which has one point per line.
   *
   * Dimension and number of points should have been
   * assigned a value before reaching this function.
   *
   * @param filename - input file
   */
  void readfile(const std::string filename) {
    std::ifstream infile;

    infile.open(filename);
    if (!infile)
      std::cout << "File " << filename << " not found!" << std::endl;

    unsigned int hRead = 0;
    float coord;
    for (unsigned int n = 0; n < N && infile; ++n) {
      for (unsigned int i = 0; i < D; ++i) {
        infile >> coord;
        try {
          points.push_back(coord);
        } catch (const std::bad_alloc&) {
          std::cerr << "Out of memory when reading data, exiting..."
                    << std::endl;
          exit(1);
        }
      }
      hRead++;
    }
    if (hRead != N)
      std::cout << "ERROR, read less than " << N << " points!!\n\n";
  }

  /** \brief Read a data set, which has 'D' and 'N' at the first line.
   *	After that, every line contains a point.
   *
   * @param filename - input file
   */
  void readfile_N_D(const std::string filename) {
    std::ifstream infile;

    infile.open(filename);
    if (!infile)
      std::cout << "File " << filename << " not found!" << std::endl;

    // read first line of file: D N
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

    unsigned int hRead = 0;
    float coord;
    for (unsigned int n = 0; n < N && infile; ++n) {
      for (unsigned int i = 0; i < D; ++i) {
        infile >> coord;
        try {
          points.push_back(coord);
        } catch (const std::bad_alloc&) {
          std::cerr << "Out of memory when reading data, exiting..."
                    << std::endl;
          exit(1);
        }
      }
      hRead++;
    }
    if (hRead != N)
      std::cout << "ERROR, read less than " << N << " points!!\n\n";
  }

#ifdef RKD_PAR
  /**
   * \brief Populate a part of the vector of roots. Used when contructing the trees in parallel.
   *
   * @param roots      	    - vector of roots
   * @param from     	 	    - starting index for populating
   * @param to         	    - ending index for populating
   * @param N         	    - size of data set
   * @param N         	    - dimension of data set
   * @param split_dims      - dimensions for splitting
   * @param t               - number of highest variances
   * to take into account for the splitting.
   * @param shuffle_enable	- flag for random shuffling
   */
  void populate(std::vector<RKD < Division_Euclidean_space<T> > >* roots, const int from,
      const int to, const size_t N, const size_t D, size_t* split_dims,
      int t, bool shuffle_enable, Rotation rotate_option) {
    //std::cout << "from = " << from << ", to = " << to << std::endl;
    try {
      if(shuffle_enable)
      std::srand(std::time(0));
      for(int i = from; i < to; ++i) {
        std::vector<size_t> indices(N);
        for (size_t i = 0; i < N; ++i) {
          indices[i] = i;
        }
        if(shuffle_enable)
          std::random_shuffle(indices.begin(), indices.end());
        (*roots)[i].build(0, indices.begin(), indices.end(), split_dims, t, D);
        if(rotate_option == Householder)
          (*roots)[i].discard_rot_p(&points);
      }
    }
    catch (const std::bad_alloc&) {
      std::cout << "[exception caught when constructing tree]\n";
      return;
    }
  }
#endif

  /**
   * \brief Search for k nearest neighbors.
   *
   * This function implements the Fast searching,
   * subsection 5.3 in the 'RKD_forest.pdf'. Notice
   * that the approximation parameters may affect
   * each other. For example, `epsilon` equal to zero,
   * can be affected by a very small value for
   * `max_leaf_check`. Moreover, the approximation
   * parameters can affect `k` as well. If for example
   * `max_leaf_check` is set to 1 and every leaf contains
   * 1 point, then we will find one NN, even if `k` >= 2.
   * In such a situation, the index of neighbors not found
   * is set to -1.
   *
   * @param q                   - query point
   * @param res                 - vector of k nearest neighbors,
   * with the first field to be the distance from
   * 'q', while the second field is the index of
   * the NN.
   * @param max_leaf_check      - max number of leaves
   * to check while searching
   * @param squared_coords      - squared coordinates of
   * points
   * @param q_squared_coords    - squared coordinates of
   * queries
   * @param q_index             - index of query
   * @param sorted_results      - if true, results are
   * going to be sorted, based on their distance
   * from the query point.
   * @param k                   - number of NN we search for
   * @param epsilon             - approximation parameter
   *
   */
  void search_nn_prune(const std::vector<T>& q,
                       std::vector<std::pair<float, int> >& res,
                       int max_leaf_check, const std::vector<T>& squared_coords,
                       const std::vector<T>& q_squared_coords, int q_index,
                       bool sorted_results = false, int k = 1, double epsilon =
                           0.0) {

    float mul_factor = 1 / (1 + epsilon);
    size_t max_count = max_leaf_check * mul_factor;
    res.resize(k, { std::numeric_limits<float>::max(), -1 });
    std::make_heap(res.begin(), res.end());
    Min_heap branch;
    size_t vis[N];
    memset(vis, 0, sizeof(vis));
    for (size_t i = 0; i < roots.size(); ++i) {
      roots[i].search_nn(q, mul_factor, max_leaf_check, res, branch, 0, 0, i,
                         vis, squared_coords, q_squared_coords, q_index);
    }
    //std::cout << "Heap.size " << branch.size() << " \n";
    size_t c = 0;
    size_t node_i, tree_i;
    float new_mindist;

    while (branch.size() && c++ < max_count) {
      std::tie(new_mindist, node_i, tree_i) = branch.top();
      branch.pop();
      //std::cout << new_mindist << std::endl;
      roots[tree_i].search_nn(q, mul_factor, max_leaf_check, res, branch,
                              new_mindist, node_i, tree_i, vis, squared_coords,
                              q_squared_coords, q_index);
    }
    //std::cout << "Checked " << c << " branches\n";
    //std::cout << "End Heap.size " << branch.size() << " \n";
    if (sorted_results)
      std::sort_heap(res.begin(), res.end());
  }

 public:
  /**
   * \brief Constructor that will automatically compute the optimal
   * parameters, build the forest, perform the search and fill in the
   * result vector. 'N' and 'D' parameters are going to be updated if
   * needed (pass them with initial values of zero to the constructor).
   *
   * @param N           - size of data set
   * @param D           - dimension of data set
   * @param datafile    - file that contains the data set
   * @param Q          	- number of queries
   * @param queryfile   - file that contains the queries
   * @param k          	- number of nearest neighbours we are searching per query
   * @param epsilon   	- factor of accuracy. Set to 0 for maximum accuracy
   * @param results			- 2D vector that will host the indices of the nearest neighbors.
   * The vector shall be resized automatically to 'Q' x 'k'.
   * @param parameters  - struct of parameters. By default, parameters are auto configured.
   * @param file_option - how to parse the `datafile`. Default value is 0.
   */
  Auto_random_kd_forest(
      size_t& N, size_t& D, const std::string& datafile, const size_t Q,
      const std::string& queryfile, const int k, const double epsilon,
      std::vector<std::vector<std::pair<float, int> > >& results,
      Params* parameters = 0,
      const int file_option = 0)
      : N(N),
        D(D) {
    // read the data set
    if (file_option == 0) {
      readfile(datafile);
    } else if (file_option == 1) {
      readfile_N_D(datafile);
    }
    assert(Auto_random_kd_forest::N > 0);
    assert(Auto_random_kd_forest::D > 0);
    N = Auto_random_kd_forest::N;
    D = Auto_random_kd_forest::D;

    float var[D];
    int points_per_leaf, trees_no, t, max_leaf_check;
    size_t sample_size;
    Rotation rotate_option;
    bool shuffle_enable;
    if(parameters) {
      assign(*parameters, points_per_leaf, trees_no, t, max_leaf_check,
             rotate_option, shuffle_enable, sample_size);
      compute_variances<T>(points, N, D, var, parameters->sample_size);
    } else {
      // variance is computed inside that function
      auto_config_params(epsilon, var, points_per_leaf, trees_no,
                         t, max_leaf_check, rotate_option, shuffle_enable);
    }

    // build the forest
    size_t split_dims[t];
    kthLargest(var, D, t, split_dims);
    /***** Compute number of leaves, nodes and maximum number of points per leaf.*****/
    size_t leaves = (N / points_per_leaf) + !!(N % points_per_leaf);  // Add one if nonzero remainder
    // If 'leaves' is a power of two, (leaves & (leaves - 1)) is zero.
    if (leaves & (leaves - 1)) {
      // Not a power of two. Round up.
      size_t old;

      // Fill all low bits
      do {
        old = leaves;
        leaves |= leaves >> 1;
      } while (old != leaves);

      // Add one
      leaves++;
    }

    size_t nmax = N / leaves + !!(N % leaves);  // +1  if remainder nonzero

    // There are one less nodes than (the power of two) number of leaves
    size_t nodes = leaves - 1;
    /************************** Done *****************************************/
#ifdef RKD_PAR
    if(rotate_option == Householder) {
      std::vector<T> rot_points[trees_no];
      for (int i = 0; i < trees_no; ++i) {
        // 'rot_points' has the scope of the for loop, watch out!
        // vector used by rotation process
        std::vector<float> v;
        rotate(N, D, rot_points[i], points, v);
        roots.push_back(RKD< Division_Euclidean_space<T> >(&rot_points[i], D, nodes, leaves, nmax, &v));
      }
      int p;
#if __GNUC__ >= 5 || __GNUC_MINOR__ >= 8
      const int P = std::thread::hardware_concurrency();
      p = (trees_no < P) ? trees_no : P;
#else
      const int P = my_hardware_concurrency();
      p = (trees_no < P) ? trees_no : P;
#endif
      if(p <= 0) {
        std::cout << "ERROR: 'p' <= 0.\n";
        return;
      }
      const int chunk = trees_no / p;

      //std::cout << "chunk = " << chunk << ", trees_no = " << trees_no << ", p = " << p <<std::endl;

      std::thread threads[p];
      // spawn 'p' threads:
      threads[p - 1] = std::thread(&Auto_random_kd_forest::populate, this, &roots, (p - 1) * chunk, trees_no,
          N, D, &split_dims[0], t, shuffle_enable, rotate_option);
      for (int i = 0; i < p - 1; ++i)
        threads[i] = std::thread(&Auto_random_kd_forest::populate, this, &roots, i * chunk, (i + 1) * chunk,
          N, D, &split_dims[0], t, shuffle_enable, rotate_option);
      for (int i = 0; i < p; ++i)
        threads[i].join();
    } else {
      for (int i = 0; i < trees_no; ++i) {
        roots.push_back(RKD< Division_Euclidean_space<T> >(&points, D, nodes, leaves, nmax));
      }
      int p;
#if __GNUC__ >= 5 || __GNUC_MINOR__ >= 8
      const int P = std::thread::hardware_concurrency();
      p = (trees_no < P) ? trees_no : P;
#else
      const int P = my_hardware_concurrency();
      p = (trees_no < P) ? trees_no : P;
#endif
      if(p <= 0) {
        std::cout << "ERROR: 'p' <= 0.\n";
        return;
      }
      const int chunk = trees_no / p;

      //std::cout << "chunk = " << chunk << ", trees_no = " << trees_no << ", p = " << p <<std::endl;

      std::thread threads[p];
      // spawn 'p' threads:
      threads[p - 1] = std::thread(&Auto_random_kd_forest::populate, this, &roots, (p - 1) * chunk, trees_no,
          N, D, &split_dims[0], t, shuffle_enable, rotate_option);
      for (int i = 0; i < p - 1; ++i)
        threads[i] = std::thread(&Auto_random_kd_forest::populate, this, &roots, i * chunk, (i + 1) * chunk,
          N, D, &split_dims[0], t, shuffle_enable, rotate_option);
      for (int i = 0; i < p; ++i)
        threads[i].join();
    }
#else // serial construction of trees
    if (shuffle_enable)
      std::srand(std::time(0));
    if (rotate_option == Householder) {
      for (int i = 0; i < trees_no; ++i) {
        std::vector<T> rot_points;
        // vector used by rotation process
        std::vector<float> v;
        rotate(N, D, rot_points, points, v);
        roots.push_back(
            RKD<Division_Euclidean_space<T> >(&rot_points, D, nodes, leaves,
                                              nmax, &v));
        std::vector<size_t> indices(N);
        for (size_t i = 0; i < N; ++i) {
          indices[i] = i;
        }
        if (shuffle_enable)
          std::random_shuffle(indices.begin(), indices.end());
        roots[i].build(0, indices.begin(), indices.end(), split_dims, t, D);
        roots[i].discard_rot_p(&points);
      }
    } else {
      for (int i = 0; i < trees_no; ++i) {
        roots.push_back(
            RKD<Division_Euclidean_space<T> >(&points, D, nodes, leaves, nmax));
        std::vector<size_t> indices(N);
        for (size_t i = 0; i < N; ++i) {
          indices[i] = i;
        }
        if (shuffle_enable)
          std::random_shuffle(indices.begin(), indices.end());
        roots[i].build(0, indices.begin(), indices.end(), split_dims, t, D);
      }
    }
#endif
    std::vector<T> squared_coords;
    computeSquare(points, N, D, squared_coords);

    std::cout << "N = " << N << ", D = " << D << std::endl;
    std::cout << "points_per_leaf = " << points_per_leaf << "\n";
    std::cout << "trees_no = " << trees_no << "\n";
    std::cout << "t = " << t << "\n";
    std::cout << "rotate_option = " << rotate_option << "\n";
    std::cout << "max_leaf_check = " << max_leaf_check << "\n";
    std::cout << std::endl;

    // read the queries
    std::vector<std::vector<T> > q;
    read_points<T>(q, Q, D, queryfile.c_str());

    // perform the search
    results.resize(Q);
    std::vector<T> q_squared_coords;
    computeSquare(q, Q, D, q_squared_coords);
    for (unsigned int i = 0; i < Q; ++i) {
      search_nn_prune(q[i], results[i], max_leaf_check, squared_coords,
                      q_squared_coords, i, false, k, epsilon);
    }
  }

  /**
   * \brief Given a match file that contains the indices
   * of the exact Nearest Neighbours ('k' integers per line,
   * total 'Q' lines) and the 'results' vector the search
   * returned, this function will check to see the quality
   * of the search, i.e. how many exact neighbours were found
   * and  it will print the miss rate.
   *
   * @param matchfile   - file of the indices of the exact NN's
   * @param Q          	- number of queries
   * @param k          	- number of nearest neighbours we searched per query
   * @param results			- 2D vector that hosts the indices of the nearest neighbors.
   */
  void check_miss(
      const std::string& matchfile, const size_t Q, const int k,
      const std::vector<std::vector<std::pair<float, int> > >& results) {
    std::vector<std::vector<int> > match;
    read_indices(match, Q, k, matchfile.c_str());
    int miss = 0;
    bool found;
    for (unsigned int i = 0; i < Q; ++i) {
      found = false;
      for (std::vector<std::pair<float, int> >::const_iterator it = results[i]
          .begin(); it != results[i].end(); ++it) {
        //std::cout << match[i][0] <<", sam: " << it->second << "\n";
        if (match[i][0] == it->second) {  // || match[i][1] == it->second) {
          found = true;
        }
        //else std::cout << "query = " << i << ", " << match[i][0] <<", " << match[i][1] << " , sam: " << it->second << "\n";
      }
      if (!found) {
        ++miss;
      }
      //    std::cout << it->first << ", index = " << it->second << std::endl;
    }

    std::cout << "miss = " << miss << std::endl;
    std::cout << "Miss: " << ((miss * 100) / (double) Q) << "%\n";
  }

  /**
   * \brief Print a point.
   *
   * @param p - the point
   */
  void print(const std::vector<T>& p) {
    for (size_t i = 0; i < p.size(); ++i)
      std::cout << p[i] << "\n";
  }

  /**
   * \brief Print the trees of the forest.
   *
   * @param points - if true, the coordinates
   * of the points are going to be printed.
   * Otherwise, the indices.
   */
  void print_trees(bool points = false) {
    for (size_t i = 0; i < roots.size(); ++i) {
      roots[i].printNodes();
      if (points)
        roots[i].printLeavesPoints();
      else
        roots[i].printLeaves();
    }
  }

 private:
  void auto_config_params(const double epsilon, float* var, int& points_per_leaf, int &trees_no,
                          int& t, int& max_leaf_check, Rotation& rotate_option,
                          bool& shuffle_enable) {
    compute_variances<T>(points, N, D, var, N);
    /* In probability theory and statistics, variance measures
     how far a set of numbers is spread out. A variance of
     zero indicates that all the values are identical. So, we
     check the first 5 variances and take their mean. If it is
     less than 1000, we choose different parameters than in the
     case that the mean is greater than 1000.
     */
    float mean_var = .0;
    for (unsigned int i = 0; i < 5 && i < D; ++i)
      mean_var += var[i];
    mean_var /= 5.0;
    //std::cout << "mean_var = " << mean_var << std::endl;

    // this should be updated with the automatic routine,
    // now hard-coded for simplicity
    if (D < 120) {
      if (N < 50000) {
        // should be updated
        if (mean_var <= 1000.0) {  //sphere-like
          points_per_leaf = 128;
          trees_no = 128;
          t = 64;
          rotate_option = No;
          shuffle_enable = false;
          max_leaf_check = 128;
          if (epsilon != 0) {
            if (epsilon <= 0.1) {
              max_leaf_check /= 4;
            } else if (epsilon <= 0.5) {
              max_leaf_check /= 8;
            } else if (epsilon > 0.5) {
              max_leaf_check /= 128;
            }
          }
        } else {  // mean of variance greater than 1000
          //3228.14 (N = 10^4)
          max_leaf_check = 1;
          trees_no = 64;
          t = 3;  //128
          rotate_option = No;
          shuffle_enable = false;
          // needs check
          points_per_leaf = 64;
          if (epsilon != 0) {
            if (epsilon <= 0.1) {
              points_per_leaf /= 2;
            } else if (epsilon <= 0.5) {
              points_per_leaf /= 4;
            } else if (epsilon > 0.5)  {
              points_per_leaf /= 8;
            }
          }
        }
      } else {  // N >= 50000 && D < 120
        points_per_leaf = 32;
        trees_no = 1024;
        t = 5; //16
        rotate_option = No;
        shuffle_enable = false;  //true
        max_leaf_check = 32;
        if (epsilon != 0) {
          if (epsilon <= 0.1) {
            max_leaf_check /= 2;
          } else if (epsilon <= 0.5) {
            max_leaf_check /= 8;
          } else if (epsilon > 0.5)  {
            max_leaf_check /= 32;
          }
        }
      }
    } else if (D < 900) {  // D >= 120 && D < 900
      // should be updated with mean variance
      // in order to pick parameters
      if (N < 100000) {
        points_per_leaf = 64;
        trees_no = 11;
        t = 64;
        rotate_option = No;
        shuffle_enable = false;
        max_leaf_check = 40;
        if (epsilon != 0) {
          if (epsilon <= 0.1) {
            max_leaf_check /= 2;
          } else if (epsilon <= 0.5) {
            max_leaf_check /= 8;
          } else if (epsilon > 0.5)  {
            max_leaf_check /= 32;
          }
        }
      } else {  // N >= 100000
        points_per_leaf = 64;
        trees_no = 128;
        t = 64;
        rotate_option = No;
        shuffle_enable = true;
        max_leaf_check = 1024;
        if(epsilon != 0) {
          if (epsilon <= 0.1) {
            max_leaf_check = 1024;
          } else if (epsilon <= 0.5) {
            max_leaf_check /= 2;
          } else if (epsilon > 0.5)  {
            max_leaf_check /= 32;
          }
        }
      }
    } else if (D < 5000) {  // D > 900 && D < 5000
      // should be updated with mean variance
      // in order to pick parameters
      if (N < 100000) {
        points_per_leaf = 64;
        trees_no = 2048;
        t = 128;
        rotate_option = No;
        shuffle_enable = true;
        max_leaf_check = 128;
        if(epsilon != 0) {
          if (epsilon <= 0.1) {
            max_leaf_check /= 2;
          } else if (epsilon <= 0.5) {
            points_per_leaf /= 2;
            max_leaf_check /= 4;
          } else if (epsilon > 0.5)  {
            points_per_leaf /= 8;  //16
            max_leaf_check /= 64;
          }
        }
      } else {
        points_per_leaf = 1024;
        trees_no = 256;
        t = 128;
        rotate_option = No;
        shuffle_enable = true;
        max_leaf_check = 128;
        if(epsilon != 0) {
          if (epsilon <= 0.1) {
            max_leaf_check /= 2;
          } else if (epsilon <= 0.5) {
            points_per_leaf /= 2;
            max_leaf_check /= 4;
          } else if (epsilon > 0.5)  {
            points_per_leaf /= 4;
            max_leaf_check /= 8;
          }
        }
      }
      if (D == 1000) {
        if (N == 1000) {
          if (mean_var > 200) {
            // klein-like data set
            points_per_leaf = 8;
            trees_no = 3;
            t = 3;  //128;
            rotate_option = No;
            shuffle_enable = true;
            max_leaf_check = 1;
            // special case
            if (epsilon != 0) {
              if (epsilon <= 0.5) {
                points_per_leaf = 6;
              } else if (epsilon > 0.5)  {
                points_per_leaf = 4;
              }
            }
          } else {
            points_per_leaf = 8;
            trees_no = 64;
            t = 8;
            rotate_option = No;
            shuffle_enable = true;
            max_leaf_check = 128;
            if (epsilon != 0) {
              if (epsilon <= 0.1) {
                max_leaf_check /= 2;
              } else if (epsilon <= 0.5) {
                points_per_leaf = 6;
                max_leaf_check /= 4;
              } else if (epsilon > 0.5)  {
                points_per_leaf = 4;
                max_leaf_check /= 8;
              }
            }
          }
        }
        if (N == 10000) {
          if (mean_var > 2000) {
            // klein-like data set
            points_per_leaf = 4;
            trees_no = 256;  //3;
            t = 3;
            rotate_option = No;
            shuffle_enable = true;
            max_leaf_check = 1;
            // special case
            if (epsilon != 0) {
              if (epsilon <= 0.5) {
                points_per_leaf = 6;
              } else if (epsilon > 0.5)  {
                points_per_leaf = 4;
              }
            }
          } else {
            points_per_leaf = 16;
            trees_no = 1024;
            t = 16;
            rotate_option = No;
            shuffle_enable = true;
            max_leaf_check = 1;
            // special case
            if (epsilon != 0) {
              if (epsilon <= 0.1) {
                points_per_leaf = 4;
              } else if (epsilon <= 0.5) {
                points_per_leaf = 2;
              } else if (epsilon > 0.5)  {
                points_per_leaf = 1;
              }
            }
          }
        }
        if (N == 100000) {  //1865.42
          if (mean_var > 1500) {
            // klein-like data set
            points_per_leaf = 64;
            trees_no = 1024;  //3;
            t = 3;
            rotate_option = No;
            shuffle_enable = true;
            max_leaf_check = 128;
            if (epsilon != 0) {
              if (epsilon <= 0.1) {
                max_leaf_check /= 2;
              } else if (epsilon <= 0.5) {
                points_per_leaf /= 2;
                max_leaf_check /= 2;
              } else if (epsilon > 0.5)  {
                points_per_leaf /= 4;
                max_leaf_check /= 4;
              }
            }
          } else {  //mean_var = 0.376393
            points_per_leaf = 64;
            trees_no = 1024;
            t = 256;
            rotate_option = No;
            shuffle_enable = true;
            max_leaf_check = 128;
            if (epsilon != 0) {
              if (epsilon <= 0.1) {
                points_per_leaf /= 2;
                max_leaf_check /= 2;
              } else if (epsilon <= 0.5) {
                points_per_leaf /= 4;
                max_leaf_check /= 4;
              } else if (epsilon > 0.5)  {
                points_per_leaf /= 16;
                max_leaf_check /= 16;
              }
            }
          }
        }
      }
    } else {  // D > 5000
     // needs further runs on data sets to improve
      points_per_leaf = 4;
      trees_no = 10;
      t = 7;
      rotate_option = No;
      shuffle_enable = false;
      max_leaf_check = 2;
      if(epsilon != 0) {
        if (epsilon <= 0.5) {
          points_per_leaf /= 2;
          max_leaf_check /= 2;
        } else if (epsilon > 0.5)  {
          points_per_leaf /= 1;
          max_leaf_check /= 1;
        }
      }
    }

    // we may have a tiny dataset, so check
    if(t > (int)D) {
      t = D - 1;
    }
    if(points_per_leaf > (int)N) {
      points_per_leaf = N;
    }
  }

  /**
   * The roots of the trees the forest has.
   */
  std::vector<RKD<Division_Euclidean_space<T> > > roots;
  /**
   * number of points
   */
  size_t N;
  /**
   * dimension of points
   */
  size_t D;
  /**
   * vector of points
   * Note that indexing is of the form: [i * D + j]
   */
  std::vector<T> points;

};
#endif /*AUTO_RANDOM_KD_TREE_H*/
