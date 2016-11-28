/**
 @file Random_kd_forest.h
 */

#ifndef RANDOM_KD_TREE_H
#define RANDOM_KD_TREE_H

#include <algorithm>
#include <tuple>
#include <queue>
#include <unordered_set>

#ifdef RKD_PAR
#include <thread>         // std::thread
#if !(__GNUC__ >= 5 || __GNUC_MINOR__ >= 8)
#include "Hardware_con.h"
#endif



// DO I NEED THESE?
#include <mutex>          // std::mutex, std::lock_guard
static std::mutex mtx;
#endif

#include "Division_Euclidean_space.h"
#include "Householder.h"
#include "Tree.h"

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
 * Enum for rotation. 
 * No means no rotation.
 * Householder means apply rotation by using
 * Householder matrices.
 * Random means just do a random shuffle of the points.
 */
enum Rotation { No, Householder, Random };

/**
 * A forest of Tree's. The points of the trees,
 * lie in the DivisionSpace.
 */
template<class DivisionSpace>
class Random_kd_forest {
  typedef typename DivisionSpace::FT FT;
#ifdef RKD_WITH_BOOST
  /**
   * The Graph we are using.
   */
  typedef typename Division_Euclidean_space<FT>::Graph Graph;
  /**
   * The vertex descriptor.
   */
  typedef typename Division_Euclidean_space<FT>::vertex_t vertex_t;
#endif

 private:
#ifdef RKD_PAR
  /**
   * \brief Populate a part of the vector of roots. Used when contructing the trees in parallel.
   *
   * @param roots      - vector of roots
   * @param from     	 - starting index for populating
   * @param to         - ending index for populating
   * @param ds         - division space, containing the data set
   * @param split_dims  - dimensions for splitting
   * @param t           - number of highest variances
   * to take into account for the splitting.
   */
  void populate(std::vector<RKD <DivisionSpace> >* roots, const int from, const int to,
                        const DivisionSpace* ds, size_t* split_dims, int t) {
    //std::cout << "from = " << from << ", to = " << to << "\n";
    try {

      std::vector<size_t> indices;
      indices.reserve(ds->size());

      for (size_t i = 0; i < ds->size(); ++i) {
        indices.push_back(i);
      }
      // using a local lock_guard to lock mtx guarantees unlocking on destruction / exception:
      //std::lock_guard<std::mutex> lck (mtx);
      for(int i = from; i < to; ++i) {
        (*roots)[i].build(0, indices.begin(), indices.end(), split_dims, t, ds->dim());
      }
    }
    catch (const std::bad_alloc&) {
      std::cout << "[exception caught when constructing tree]\n";
      return;
    }
  }
#endif

 public:

  /**
   * \brief Constructor that will build the trees, indexing the points of
   * the given DivisionSpace.
   *
   * @param ds                - the DivisionSpace
   * @param sample            - sample size when computing variances
   * @param points_per_leaf   - desired number of points per leaf
   * @param trees_no          - number of trees
   * @param t                 - number of dimensions to be
   * used for building the tree. 5 seems to work well.
   * @param rotate_option   	- No, Householder or Random. Default is No.
   */
  Random_kd_forest(DivisionSpace & ds, const int sample = 100,
                 size_t points_per_leaf = 1, const int trees_no = 5, int t = 5,
                 Rotation rotate_option = No) {

    /***** Compute number of leaves, nodes and maximum number of points per leaf.*****/
    size_t leaves = (ds.size() / points_per_leaf)
        + !!(ds.size() % points_per_leaf);  // Add one if nonzero remainder
    /*
     leaves = ( forest->points.size() / points_per_leaf )
     + ( dataset->points % perleaf > 0) ? 1 : 0;
     */

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

    size_t nmax = ds.size() / leaves + !!(ds.size() % leaves);  // +1  if remainder nonzero

    // There are one less nodes than (the power of two) number of leaves
    size_t nodes = leaves - 1;

    /************************** Done *****************************************/

    //std::cout << "nodes = " << nodes << ", nmax = " << nmax
    //<< ", leaves = " << leaves << ", total points = " << points.size() << "\n";
    if (rotate_option == Householder) {
      DivisionSpace rot_dataset(ds.size(), ds.dim());
      // vector used by rotation process
      std::vector<float> v;

      rotate(ds.size(), ds.dim(), rot_dataset, ds.getPoints(), v);

      float avg[ds.dim()];
      float var[ds.dim()];

      size_t split_dims[t];

      compute_variances<FT>(t, ds.getPoints(), ds.dim(), ds.size(), avg, var,
                        split_dims, sample);

#ifdef RKD_PAR
      for (int i = 0; i < trees_no; ++i) {
              roots.push_back(
                  RKD<DivisionSpace>(&rot_dataset, nodes, leaves, nmax, &v, &ds));
      }
      /*
       *  If your vector has n elements and you have p threads,
       *  thread i writes only to elements [i n / p, (i + 1) n / p).
       *  Note that this is preferable over having thread i write to
       *  elements at index j only if j mod p = i because it leads to
       *  fewer cache invalidations.
       */
      int p;
      #if __GNUC__ >= 5 || __GNUC_MINOR__ >= 8 // 4.8 for example
        const int P = std::thread::hardware_concurrency();
      	p = (trees_no < P) ? trees_no : P;
      	//std::cout << P << " concurrent threads are supported.\n";
      #else
      	const int P = my_hardware_concurrency();
        p = (trees_no < P) ? trees_no : P;
        //std::cout << P << " concurrent threads are supported.\n";
      #endif
      if(!p) {
      	std::cout << "ERROR: 'p' is zero.\n";
      	return;
      }
      const int chunk = trees_no / p;

      std::thread threads[p];
      // spawn 'p' threads:
      threads[p - 1] = std::thread(&Random_kd_forest::populate, this, &roots, (p - 1) * chunk, trees_no,
                                   &rot_dataset, &split_dims[0], t);
      for (int i = 0; i < p - 1; ++i)
        threads[i] = std::thread(&Random_kd_forest::populate, this, &roots, i * chunk, (i + 1) * chunk,
                                 &rot_dataset, &split_dims[0], t);

      for (int i = 0; i < p; ++i)
      	threads[i].join();
#else // serial construction of trees
      for (int i = 0; i < trees_no; ++i) {
        roots.push_back(
            RKD<DivisionSpace>(&rot_dataset, nodes, leaves, nmax, &v, &ds));
        std::vector<size_t> indices;
        indices.reserve(ds.size());

        for (size_t i = 0; i < ds.size(); ++i) {
          indices.push_back(i);
        }
        roots[i].build(0, indices.begin(), indices.end(), split_dims, t, ds.dim());
      }
#endif
    } else {
      float avg[ds.dim()];
      float var[ds.dim()];

      size_t split_dims[t];

      compute_variances<FT>(t, ds.getPoints(), ds.dim(), ds.size(), avg, var,
                        split_dims, sample);

#ifdef RKD_PAR
      for (int i = 0; i < trees_no; ++i) {
              roots.push_back(
                  RKD<DivisionSpace>(&ds, nodes, leaves, nmax, nullptr, nullptr));
      }
      /*
       *  If your vector has n elements and you have p threads,
       *  thread i writes only to elements [i n / p, (i + 1) n / p).
       *  Note that this is preferable over having thread i write to
       *  elements at index j only if j mod p = i because it leads to
       *  fewer cache invalidations.
       */    
      int p;
      #if __GNUC__ >= 5 || __GNUC_MINOR__ >= 8 // 4.8 for example
        const int P = std::thread::hardware_concurrency();
      	p = (trees_no < P) ? trees_no : P;
      	//std::cout << P << " concurrent threads are supported.\n";
      #else
      	const int P = my_hardware_concurrency();
        p = (trees_no < P) ? trees_no : P;
        //std::cout << P << " concurrent threads are supported.\n";
      #endif
      if(!p) {
      	std::cout << "ERROR: 'p' is zero.\n";
      	return;
      }
      const int chunk = trees_no / p;

      std::thread threads[p];
      // spawn 'p' threads:
      threads[p - 1] = std::thread(&Random_kd_forest::populate, this, &roots, (p - 1) * chunk, trees_no,
                                   &ds, &split_dims[0], t);
      for (int i = 0; i < p - 1; ++i) 
        threads[i] = std::thread(&Random_kd_forest::populate, this, &roots, i * chunk, (i + 1) * chunk,
                                 &ds, &split_dims[0], t);

      for (int i = 0; i < p; ++i)
      	threads[i].join();
#else // serial construction of trees
      for (int i = 0; i < trees_no; ++i) {
        roots.push_back(
            RKD<DivisionSpace>(&ds, nodes, leaves, nmax, nullptr, nullptr));
        std::vector<size_t> indices;
        indices.reserve(ds.size());

        for (size_t i = 0; i < ds.size(); ++i) {
          indices.push_back(i);
        }
        if(rotate_option == Random)
        	std::random_shuffle(indices.begin(), indices.end());
        roots[i].build(0, indices.begin(), indices.end(), split_dims, t, ds.dim());
      }
#endif
    }
  }

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
   * to check while searching/
   * @param sorted_results      - if true, results are
   * going to be sorted, based on their distance
   * from the query point.
   * @param k                   - number of NN we search for
   * @param epsilon             - approximation parameter
   *
   */
  void search_nn_prune(const std::vector<FT>& q,
                       std::vector<std::pair<float, int> >& res,
                       int max_leaf_check, bool sorted_results = false, int k =
                           1,
                       float epsilon = 0) {
                       
    float mul_factor = 1 / (1 + epsilon);
    size_t max_count = max_leaf_check * mul_factor;
    res.resize(k, { std::numeric_limits<float>::max(), -1 });
    std::make_heap(res.begin(), res.end());
    Min_heap branch;
    std::unordered_set<size_t> vis;
    for (size_t i = 0; i < roots.size(); ++i) {
      roots[i].search_nn(q, mul_factor, max_leaf_check, res, branch, 0, 0, i, vis);
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
                              new_mindist, node_i, tree_i, vis);
    }
    //std::cout << "Checked " << c << " branches\n";
    //std::cout << "End Heap.size " << branch.size() << " \n";
    if (sorted_results)
      std::sort_heap(res.begin(), res.end());
  }
  
#ifdef RKD_WITH_BOOST
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
   * @param g                   - the graph
   * @param v										- vertex descriptor
   * @param res                 - vector of k nearest neighbors,
   * with the first field to be the distance from
   * 'q', while the second field is the index of
   * the NN.
   * @param max_leaf_check      - max number of leaves
   * to check while searching/
   * @param sorted_results      - if true, results are
   * going to be sorted, based on their distance
   * from the query point.
   * @param k                   - number of NN we search for
   * @param epsilon             - approximation parameter
   *
   */
  void search_nn_prune(Graph& g, vertex_t& v,
                       std::vector<std::pair<float, int> >& res,
                       int max_leaf_check, bool sorted_results = false, int k =
                           1,
                       float epsilon = 0) {

		std::vector<FT> q;
		for(unsigned int i = 0; i < g[v].v.size(); ++i)
			q.push_back(g[v].v[i]);
    float mul_factor = 1 / (1 + epsilon);
    size_t max_count = max_leaf_check * mul_factor;
    res.resize(k, { std::numeric_limits<float>::max(), -1 });
    std::make_heap(res.begin(), res.end());
    Min_heap branch;
    std::unordered_set<size_t> vis;
    for (size_t i = 0; i < roots.size(); ++i) {
      roots[i].search_nn(q, mul_factor, max_leaf_check, res, branch, 0, 0, i, vis);
    }
    size_t c = 0;
    size_t node_i, tree_i;
    float new_mindist;

    while (branch.size() && c++ < max_count) {
      std::tie(new_mindist, node_i, tree_i) = branch.top();
      branch.pop();
      roots[tree_i].search_nn(q, mul_factor, max_leaf_check, res, branch,
                              new_mindist, node_i, tree_i, vis);
    }
    if (sorted_results)
      std::sort_heap(res.begin(), res.end());
  }
#endif

  /**
   * \brief Search for k nearest neighbors, with every
   * point of the data set to play the role of the query.
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
   * @param N                   - number of points
   * @param res                 - vector of vectors
   * of k nearest neighbors, with the first field
   * to be the distance from 'q', while the second
   * field is the index of the NN. The size of the
   * vector is equal to the size of the data set.
   * @param max_leaf_check      - max number of leaves
   * to check while searching/
   * @param sorted_results      - if true, results are
   * going to be sorted, based on their distance
   * from the query point.
   * @param k                   - number of NN we search for
   * @param epsilon             - approximation parameter
   *
   */
  void search_nn_prune(
      const size_t N,
      std::vector<std::vector<std::pair<float, int> > >& res,
      int max_leaf_check, bool sorted_results, int k = 1, float epsilon = 0) {

    float mul_factor = 1 / (1 + epsilon);
    size_t max_count = max_leaf_check * mul_factor;
    const std::vector<FT>* points = roots[0].get_points();
    const size_t D = roots[0].dim();
    res.resize(N);

    for (size_t i = 0; i < N; ++i) {
      res[i].resize(k, { std::numeric_limits<float>::max(), -1 });
      std::make_heap(res[i].begin(), res[i].end());
      Min_heap branch;
      std::unordered_set<size_t> vis;
      int offset = i * D;
      std::vector<FT> q(points->begin() + offset, points->begin() + offset + D);
      for (size_t j = 0; j < roots.size(); ++j) {
        roots[j].search_nn(q, mul_factor, max_leaf_check, res[i], branch, 0, 0,
                           j, vis);
      }
      size_t c = 0;
      size_t node_i, tree_i;
      float new_mindist;

      while (branch.size() && c++ < max_count) {
        std::tie(new_mindist, node_i, tree_i) = branch.top();
        branch.pop();
        roots[tree_i].search_nn(q, mul_factor, max_leaf_check, res[i], branch,
                                new_mindist, node_i, tree_i, vis);
      }
      if (sorted_results)
        std::sort_heap(res[i].begin(), res[i].end());
    }
  }

  /**
   * \brief Search for k **exact** nearest neighbors.
   *
   * This function implements the Exact searching,
   * subsection 5.1 in the 'RKD_forest.pdf'. Notice
   * that the approximation parameters may affect
   * each other. For example, `epsilon` equal to zero,
   * can be affected by a very small value for
   * `max_leaf_check`. Moreover, the approximation
   * parameters can affect `k` as well. If for example
   * `max_leaf_check` is set to 1 and every leaf contains
   * 1 point, then we will find one NN, even if `k` >= 2.
   * In such a situation, the index of neighbors not found
   * is set to -1.
   * For **exact** searching,
   * use 'epsilon' equal to 0 and a big value
   * for `max_leaf_check` and **one tree**.
   *
   * @param q                   - query point
   * @param res                 - vector of k nearest neighbors,
   * with the first field to be the distance from
   * 'q', while the second field is the index of
   * the NN.
   * @param max_leaf_check      - max number of leaves
   * to check while searching/
   * @param k                   - number of NN we search for
   * @param epsilon             - approximation parameter
   *
   */
  void search_nn(const std::vector<FT>& q,
                 std::vector<std::pair<float, int> >& res,
                 int max_leaf_check, unsigned int k = 1, float epsilon = 0) {
    float mul_factor = 1 / (1 + epsilon);
    std::vector<std::vector<std::pair<float, int> > > v;
    v.resize(roots.size());
    std::unordered_set<size_t> visited;
    for (size_t i = 0; i < roots.size(); ++i) {
      v[i] = std::vector<std::pair<float, int> > { k, { std::numeric_limits<
          FT>::max(), -1 } };
      std::make_heap(v[i].begin(), v[i].end());
      roots[i].search_nn(q, mul_factor, max_leaf_check, v[i], visited);
      std::sort_heap(v[i].begin(), v[i].end());
    }

    std::vector<std::pair<float, int> > con;
    con.reserve(roots.size() * k);
    for (size_t i = 0; i < roots.size(); ++i)
      con.insert(con.end(), v[i].begin(), v[i].end());
    std::sort(con.begin(), con.end(), lib_compar);

    int i = 0;
    for (int j = 0; j < k; ++i) {
      if ((std::find(res.begin(), res.end(), con[i])) == res.end()) {
        res.push_back(con[i]);
        ++j;
      }
    }
  }
  
#ifdef RKD_WITH_BOOST
  /**
   * \brief Search for k **exact** nearest neighbors.
   *
   * This function implements the Exact searching,
   * subsection 5.1 in the 'RKD_forest.pdf'. Notice
   * that the approximation parameters may affect
   * each other. For example, `epsilon` equal to zero,
   * can be affected by a very small value for
   * `max_leaf_check`. Moreover, the approximation
   * parameters can affect `k` as well. If for example
   * `max_leaf_check` is set to 1 and every leaf contains
   * 1 point, then we will find one NN, even if `k` >= 2.
   * In such a situation, the index of neighbors not found
   * is set to -1.
   * For **exact** searching,
   * use 'epsilon' equal to 0 and a big value
   * for `max_leaf_check` and **one tree**.
   *
   * @param g                   - the graph
   * @param v_d									- vertex descriptor
   * @param res                 - vector of k nearest neighbors,
   * with the first field to be the distance from
   * 'q', while the second field is the index of
   * the NN.
   * @param max_leaf_check      - max number of leaves
   * to check while searching/
   * @param k                   - number of NN we search for
   * @param epsilon             - approximation parameter
   *
   */
  void search_nn(Graph& g, vertex_t& v_d,
                 std::vector<std::pair<float, int> >& res,
                 int max_leaf_check, int k = 1, float epsilon = 0) {
    std::vector<FT> q;
		for(int i = 0; i < g[v_d].v.size(); ++i)
			q.push_back(g[v_d].v[i]);
    float mul_factor = 1 / (1 + epsilon);
    std::vector<std::vector<std::pair<float, int> > > v;
    v.resize(roots.size());
    std::unordered_set<size_t> visited;
    for (size_t i = 0; i < roots.size(); ++i) {
      v[i] = std::vector<std::pair<float, int> > { k, { std::numeric_limits<
          FT>::max(), -1 } };
      std::make_heap(v[i].begin(), v[i].end());
      roots[i].search_nn(q, mul_factor, max_leaf_check, v[i], visited);
      std::sort_heap(v[i].begin(), v[i].end());
    }

    std::vector<std::pair<float, int> > con;
    con.reserve(roots.size() * k);
    for (size_t i = 0; i < roots.size(); ++i)
      con.insert(con.end(), v[i].begin(), v[i].end());
    std::sort(con.begin(), con.end(), lib_compar);

    int i = 0;
    for (int j = 0; j < k; ++i) {
      if ((std::find(res.begin(), res.end(), con[i])) == res.end()) {
        res.push_back(con[i]);
        ++j;
      }
    }
  }
#endif

  /**
   * \brief Search for k nearest neighbors,
   * based on a given radius.
   *
   * This function implements the Exact radius searching,
   * subsection 5.2 in the 'RKD_forest.pdf'. Notice
   * that the approximation parameters may affect
   * each other. For example, `epsilon` equal to zero,
   * can be affected by a very small value for
   * `max_leaf_check`.
   *
   * @param q                   - query point
   * @param res                 - vector of k nearest neighbors,
   * with the first field to be the index of NN,
   * while the second field is the distance from
   * 'q'.
   * @param r                   - radius
   * @param epsilon             - approximation parameter
   * @param max_leaf_check      - max number of leaves
   * to check while searching/
   *
   * @return                    - 0 if no results are found, else
   * number of results
   */
  size_t search_nn_radius(std::vector<FT>& q,
                        std::vector<std::pair<int, float> >& res, float r,
                        float epsilon = 0, int max_leaf_check = 50) {

    float mul_factor = 1 / (1 + epsilon);
    std::vector<std::vector<std::pair<int, float> > > v;
    v.resize(roots.size());
    std::unordered_set<size_t> visited;
    for (size_t i = 0; i < roots.size(); ++i) {
      roots[i].search_nn_radius(q, mul_factor, max_leaf_check, r, v[i], visited);
    }

    res.insert(res.end(), v[0].begin(), v[0].end());
    for (size_t i = 1; i < roots.size(); ++i) {
      for (size_t j = 0; j < v[i].size(); ++j) {
        if ((std::find(res.begin(), res.end(), v[i][j])) == res.end()) {
          res.push_back(v[i][j]);
        }
      }
    }
    return res.size();
  }
  
#ifdef RKD_WITH_BOOST
    /**
   * \brief Search for k nearest neighbors,
   * based on a given radius.
   *
   * This function implements the Exact radius searching,
   * subsection 5.2 in the 'RKD_forest.pdf'. Notice
   * that the approximation parameters may affect
   * each other. For example, `epsilon` equal to zero,
   * can be affected by a very small value for
   * `max_leaf_check`.
   *
   * @param g                   - the graph
   * @param v_d									- vertex descriptor
   * @param res                 - vector of k nearest neighbors,
   * with the first field to be the index of NN,
   * while the second field is the distance from
   * 'q'.
   * @param r                   - radius
   * @param epsilon             - approximation parameter
   * @param max_leaf_check      - max number of leaves
   * to check while searching/
   *
   * @return                    - 0 if no results are found, else
   * number of results
   */
  size_t search_nn_radius(Graph& g, vertex_t& v_d,
                        std::vector<std::pair<int, float> >& res, float r,
                        float epsilon = 0, int max_leaf_check = 50) {

		std::vector<FT> q;
		for(int i = 0; i < g[v_d].v.size(); ++i)
			q.push_back(g[v_d].v[i]);
    float mul_factor = 1 / (1 + epsilon);
    std::vector<std::vector<std::pair<int, float> > > v;
    v.resize(roots.size());
    std::unordered_set<size_t> visited;
    for (size_t i = 0; i < roots.size(); ++i) {
      roots[i].search_nn_radius(q, mul_factor, max_leaf_check, r, v[i], visited);
    }

    res.insert(res.end(), v[0].begin(), v[0].end());
    for (size_t i = 1; i < roots.size(); ++i) {
      for (size_t j = 0; j < v[i].size(); ++j) {
        if ((std::find(res.begin(), res.end(), v[i][j])) == res.end()) {
          res.push_back(v[i][j]);
        }
      }
    }
    return res.size();
  }
#endif
  /**
   * \brief Print a point.
   *
   * @param p - the point
   */
  void print(const std::vector<FT>& p) {
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
  /**
   * The roots of the trees the forest has.
   */
  std::vector<RKD <DivisionSpace> > roots;

};
#endif /*RANDOM_KD_TREE_H*/
