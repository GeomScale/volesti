/**
 @file Tree.h
 */

#ifndef TREE_H
#define TREE_H

#include <limits>
#include "Mean_variance.h"
#include "Householder.h"
#include "Find_diameter.h"

/**
 * Min_heap is actually a std::priority_queue,
 * with std::greater as a parameter.
 */
typedef std::priority_queue<std::tuple<float, int, int>,
    std::vector<std::tuple<float, int, int> >,
    std::greater<std::tuple<float, int, int> > > Min_heap;

/**
 * \brief RKD's nodes and leaves are
 * stored in arrays. Building and
 * searching for k-NN supported.
 *
 */
template<class DivisionSpace>
class RKD {
 public:

  /**
   * The data type we are using.
   */
  typedef typename DivisionSpace::FT FT;

  /**
   * \brief Constructor of the tree.
   *
   * @param ds          - the space where the points
   * lie into
   * @param nodes_size  - number of nodes
   * @param leaves_size - number of leaves
   * @param nmax        - maximum number of points per leaf
   * @param rot_v       - rotation vector
   * @param ds_original          - the space where the original
   * points lie into
   */
  RKD(const DivisionSpace* ds, size_t nodes_size,
      size_t leaves_size, size_t nmax, std::vector<float>* rot_v = 0,
      const DivisionSpace* ds_original = 0)
      : points(&(ds->getPoints())),
        D(ds->dim()),
        l(leaves_size),
        nmax(nmax),
        n(nodes_size) {

    // might be buggy in the rotation case

    leaves.resize(leaves_size);
    nodes.resize(nodes_size);

    if (rot_v)
      rot = *rot_v;

    if (ds_original)
      points = &(ds_original->getPoints());    
    //printNodes();
    //printLeaves();
  }
  
  /**
   * \brief Constructor of the tree.
   *
   * @param p          	- the points
   * @param D						- dimension of points
   * @param nodes_size  - number of nodes
   * @param leaves_size - number of leaves
   * @param nmax        - maximum number of points per leaf
   * @param rot_v				-	rotation vector
   */
  RKD(const std::vector<FT>* p, const size_t D, size_t nodes_size,
      size_t leaves_size, size_t nmax, std::vector<float>* rot_v = 0)
      : points(p),
        D(D),
        l(leaves_size),
        nmax(nmax),
        n(nodes_size) {

    leaves.resize(leaves_size);
    nodes.resize(nodes_size);

    if (rot_v)
      rot = *rot_v;

    //printNodes();
    //printLeaves();
  }

  /**
   * \brief Recursively build the tree.
   *
   * @param index       - index of node/leaf to be created
   * @param begin       - iterator pointing to start of the data
   * we process
   * @param end         - iterator pointing to end of the data
   * we process
   * @param split_dims  - dimensions for splitting
   * @param t           - number of highest variances
   * to take into account for the splitting.
   * @param D           - dimension of the points
   */
  void build(size_t index, std::vector<size_t>::iterator begin,
             std::vector<size_t>::iterator end, size_t* split_dims,
             int t, int D) {
    if (index >= n) {  // create leaf
      size_t leaf_index = index - n;
      for (std::vector<size_t>::iterator it = begin; it < end; ++it)
        leaves[leaf_index].push_back(*it);
      return;
    } else {  // create node
      int split_dim = random(t - 1);

//            std::cout << "before median init\n";
//            for (std::vector<size_t>::iterator it = begin; it < end; ++it)
//              std::cout << *it << " ";
//            std::cout << "\nbefore median done\n" << std::endl;

      // median holds the index of the point that is the median in
      // 'split_dims[split_dim]' dimension
      std::vector<size_t>::iterator median = find_median(
          *points, begin, end, ((end - begin) / 2), split_dims[split_dim],
          D);

//            std::cout << "median = " << *median << "\n";
//            for (std::vector<size_t>::iterator it = begin; it < end; ++it)
//              std::cout << *it << " ";
//            std::cout << std::endl;
//            std::cout << "split dim = " << split_dim << std::endl;
//            std::cout << "split_dims[split_dim] = " << split_dims[split_dim]
//                                                                  << std::endl;
//            std::cout << "cut value = " <<
//                (*points)[*median*D + split_dims[split_dim]]
//                                                     << std::endl << std::endl;

      //float diameter = find_diameter_appr(*points, begin, end, D);
      //float diameter = find_diameter_exact(*points, begin, end, D);
	    //std::cout << "DIAMETER_appr = " << diameter << "\n";
      //float bound = (3 * diameter)/sqrt(D);
      //float delta = random<float>(-bound, bound);
      //std::cout << "delta = " << delta << ", bound = " << bound << "\n";

      // Note that indexing is of the form: [i * D + j]
      nodes[index] = {split_dims[split_dim],
        (*points)[*median*D + split_dims[split_dim]] /*+delta*/};

      build(2 * index + 1, begin, begin + (median - begin), split_dims, t, D);
      build(2 * index + 2, begin + (median - begin), end, split_dims, t, D);
    }
  }

  /**
   * \brief Must be called after the building
   * process, in order to replace data member
   * 'points' with the original points and
   * discard the rotated points.
   */
  void discard_rot_p(const std::vector<FT>* p) {
    points = p;
  }

  /**
   * \brief Print nodes of tree.
   */
  void printNodes() {
    for (size_t i = 0; i < nodes.size(); ++i)
      std::cout << nodes[i].first << " " << nodes[i].second << "\n";
  }

  /**
   * \brief Print leaves of tree.
   */
  void printLeaves() {
    for (size_t i = 0; i < leaves.size(); ++i) {
      std::cout << "\n" << i << "-th leaf: ";
      for (size_t j = 0; j < nmax; ++j) {
        std::cout << leaves[i][j] << " ";
      }
    }
    std::cout << "\nDone printing leaves\n";
  }

  /**
   * \brief Print the coordinates
   * of the leaves of tree.
   */
  void printLeavesPoints() {
    for (size_t i = 0; i < leaves.size(); ++i) {
      std::cout << "\n" << i << "-th leaf: ";
      for (size_t j = 0; j < nmax; ++j) {
        std::cout << "\n";
        for (size_t k = 0; k < D; ++k) {
          std::cout << (*points)[leaves[i][j] * D + k] << " ";
        }
      }
    }
    std::cout << "\nDone printing leaves\n";
  }

  /**
   * \brief Recursively search for nearest neighbors.
   *
   * @param index           - index of node/leaf we are accessing
   * @param q_rot           - rotated query
   * @param q               - original query
   * @param mul_factor      - approximation factor
   * @param res             - vector of results (NNs)
   * @param c               - number of leaves visited
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param vis             - points already visited
   */
  void nn(size_t index, const std::vector<FT>& q_rot, const std::vector<FT>& q,
          const float& mul_factor, std::vector <std::pair<float,
          int> >& res, int& c, const int max_leaf_check,
          std::unordered_set<size_t>& vis) {
    if(c == max_leaf_check)
      return;
    if (index >= n) {
      ++c;
      // leaf
      size_t point_index, leaf_index = index - n;
      float dist;
      for (size_t j = 0; j < leaves[leaf_index].size(); ++j) {
        // leaves[i][j] gives you the index of the j-th point,
        // located in the i-th leaf.
        // (*points)[leaves[leaf_index][j]] is where the first
        // coordinate of the point described above lies.
        point_index = leaves[leaf_index][j] * D;
        if(vis.find(point_index) != vis.end()) {
          continue;
        }
        vis.insert(point_index);
        dist = Eucl_distance(*points, point_index, point_index + D, q);
        dist *= mul_factor;

        if (dist < res.front().first) {
          std::pop_heap (res.begin(), res.end());
          res.pop_back();
          res.push_back({dist, leaves[leaf_index][j]});
          std::push_heap (res.begin(), res.end());
        }
      }
    } else {
      // node
      if (q_rot[nodes[index].first] < nodes[index].second) {
        // go to left child
        nn(2 * index + 1, q_rot, q, mul_factor, res, c, max_leaf_check, vis);
        //std::cout << "f " << q_rot[nodes[index].first] << " "
        //    << best << " " << nodes[index].second << "\n";
        if (q_rot[nodes[index].first] + res.front().first >= nodes[index].second) {
          nn(2 * index + 2, q_rot, q, mul_factor, res, c, max_leaf_check, vis);
        }
      } else {
        // go to right child
        nn(2 * index + 2, q_rot, q, mul_factor, res, c, max_leaf_check, vis);
        //std::cout << "s " << q_rot[nodes[index].first] << " "
        //    << best << " " << nodes[index].second << "\n";
        if (q_rot[nodes[index].first] - res.front().first < nodes[index].second) {
          nn(2 * index + 1, q_rot, q, mul_factor, res, c, max_leaf_check, vis);
        }
      }
    }
  }

  /**
   * \brief Recursively search for nearest neighbors.
   *
   * @param index           - index of node/leaf we are accessing
   * @param q               - original query
   * @param mul_factor      - approximation factor
   * @param res             - vector of results (NNs)
   * @param c               - number of leaves visited
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param vis             - points already visited
   */
  void nn(size_t index, const std::vector<FT>& q, const float& mul_factor,
          std::vector <std::pair<float, int> >& res, int& c,
          const int max_leaf_check, std::unordered_set<size_t>& vis) {
    if(c == max_leaf_check)
      return;
    if (index >= n) {
      ++c;
      // leaf
      size_t point_index, leaf_index = index - n;
      float dist;
      for (size_t j = 0; j < leaves[leaf_index].size(); ++j) {
        // leaves[i][j] gives you the index of the j-th point,
        // located in the i-th leaf.
        // (*points)[leaves[leaf_index][j]] is where the first
        // coordinate of the point described above lies.
        point_index = leaves[leaf_index][j] * D;
        if(vis.find(point_index) != vis.end()) {
          continue;
        }
        vis.insert(point_index);
        dist = Eucl_distance(*points, point_index, point_index + D, q);
        dist *= mul_factor;

        if (dist < res.front().first) {
          std::pop_heap (res.begin(), res.end());
          res.pop_back();
          res.push_back({dist, leaves[leaf_index][j]});
          std::push_heap (res.begin(), res.end());
        }
      }
    } else {
      // node
      if (q[nodes[index].first] < nodes[index].second) {
        // go to left child
        nn(2 * index + 1, q, mul_factor, res, c, max_leaf_check, vis);
        if (q[nodes[index].first] + res.front().first >= nodes[index].second) {
          nn(2 * index + 2, q, mul_factor, res, c, max_leaf_check, vis);
        }
      } else {
        // go to right child
        nn(2 * index + 2, q, mul_factor, res, c, max_leaf_check, vis);
        if (q[nodes[index].first] - res.front().first < nodes[index].second) {
          nn(2 * index + 1, q, mul_factor, res, c, max_leaf_check, vis);
        }
      }
    }
  }

  /**
   * \brief Search for nearest neighbors.
   *
   * @param q               - original query
   * @param mul_factor      - approximation factor
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param res             - vector of results (NNs)
   * @param visited         - visited points
   */
  void search_nn(const std::vector<FT>& q, float& mul_factor, int max_leaf_check,
                 std::vector <std::pair<float, int> >& res,
                 std::unordered_set<size_t>& visited) {
    int c = 0; // number of leaves visited
    if (rot.size()) {
      std::vector<FT> q_rot;
      householder_transform(1, D, q_rot, q, rot);
      nn(0, q_rot, q, mul_factor, res, c, max_leaf_check, visited);
    } else {
      nn(0, q, mul_factor, res, c, max_leaf_check, visited);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  /**
   * \brief Recursively search for nearest neighbors.
   *
   * @param index           - index of node/leaf we are accessing
   * @param q_rot           - rotated query
   * @param q               - original query
   * @param mul_factor      - approximation factor
   * @param res             - vector of results (NNs)
   * @param c               - number of leaves visited
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param vis             - points already visited
   * @param branch          - minimum heap of unvisited branches
   * @param mindist         - current minimum distance
   * @param tree_i          - tree index
   * @param squared_coords   - squared coordinates of the points
   * @param q_squared_coords - squared coordinates of the queries
   * @param q_index          - index of query
   */
  void nn(size_t index, const std::vector<FT>& q_rot, const std::vector<FT>& q,
          const float& mul_factor, std::vector <std::pair<float,
          int> >& res, int& c, const int max_leaf_check,
          size_t* vis,
          Min_heap& branch,
          const float& mindist, const size_t& tree_i,
          const std::vector<FT>& squared_coords,
          const std::vector<FT>& q_squared_coords, int q_index) {
    if(c == max_leaf_check || res.front().first < mindist)
      return;
    if (index >= n) {
      ++c;
      // leaf
      size_t point_index, leaf_index = index - n;
      float dist;
      for (size_t j = 0; j < leaves[leaf_index].size(); ++j) {
        // leaves[i][j] gives you the index of the j-th point,
        // located in the i-th leaf.
        // (*points)[leaves[leaf_index][j]] is where the first
        // coordinate of the point described above lies.
        point_index = leaves[leaf_index][j] * D;
        if(vis[leaves[leaf_index][j]])
          continue;
        vis[leaves[leaf_index][j]] = 1;
        dist = Eucl_distance(*points, point_index, point_index + D, q);
        dist *= mul_factor;

        if (dist < res.front().first) {
          std::pop_heap (res.begin(), res.end());
          res.pop_back();
          res.push_back({dist, leaves[leaf_index][j]});
          std::push_heap (res.begin(), res.end());
        }
      }
    } else {
      // node
      FT val_q = q_rot[nodes[index].first];
      FT val_n = nodes[index].second;

      /* Which child branch should be taken first? */
      FT diff = val_q - val_n;
      size_t best_child_i, other_child_i;
      if(diff < 0) {
        best_child_i = 2 * index + 1;
        other_child_i = 2 * index + 2;
      }
      else {
        best_child_i = 2 * index + 2;
        other_child_i = 2 * index + 1;
      }

      /*
       *  Partial euclidean distance, using just one dimension. This is used by the
       *  kd-tree when computing partial distances while traversing the tree.
       *
       *  Square root is omitted for efficiency.
       */
      float new_dist = mindist + ((val_q - val_n) * (val_q - val_n));

      // collect branch not taken for future use
      if(new_dist * mul_factor < res.front().first) {
          branch.push(std::make_tuple(new_dist, other_child_i, tree_i));
      }

      nn(best_child_i, q_rot, q, mul_factor, res, c, max_leaf_check, vis,
         branch, mindist, tree_i, squared_coords, q_squared_coords, q_index);
    }
  }

  /**
   * \brief Recursively search for nearest neighbors.
   *
   * @param index            - index of node/leaf we are accessing
   * @param q                - original query
   * @param mul_factor       - approximation factor
   * @param res              - vector of results (NNs)
   * @param c                - number of leaves visited
   * @param max_leaf_check   - max number of
   * leaves to check
   * @param vis              - points already visited
   * @param branch           - minimum heap of unvisited branches
   * @param mindist          - current minimum distance
   * @param tree_i           - tree index
   * @param squared_coords   - squared coordinates of the points
   * @param q_squared_coords - squared coordinates of the queries
   * @param q_index          - index of query
   */
  void nn(size_t index, const std::vector<FT>& q, const float& mul_factor,
          std::vector <std::pair<float, int> >& res, int& c,
          const int max_leaf_check, size_t* vis,
          Min_heap& branch,
          const float& mindist, const size_t& tree_i,
          const std::vector<FT>& squared_coords,
          const std::vector<FT>& q_squared_coords, int q_index) {
    //edw benw
    if(c == max_leaf_check || res.front().first < mindist)
      return;
    if (index >= n) {
      ++c;
      // leaf
      size_t point_index, leaf_index = index - n;
      float dist;
      for (size_t j = 0; j < leaves[leaf_index].size(); ++j) {
        // leaves[i][j] gives you the index of the j-th point,
        // located in the i-th leaf.
        // (*points)[leaves[leaf_index][j]] is where the first
        // coordinate of the point described above lies.
        point_index = leaves[leaf_index][j] * D;
        //std::cout << "point_index = " << point_index/D << "\n";
        //if(point_index/D == 4096 ||point_index/D == 3138)
        //	std::cout <<"fouuuuuuuuuuuuuuuuuuuuuuuuuuuuund\n";
        if(vis[leaves[leaf_index][j]])
          continue;
        vis[leaves[leaf_index][j]] = 1;

        //dist = squared_Eucl_distance(*points, point_index, point_index + D, q);
        //std::cout << "Normal: = " << dist << std::endl;
        dist = squared_Eucl_distance(*points, point_index, point_index + D,
                                     squared_coords[point_index/D], q,
                                     q_squared_coords[q_index]);
        //std::cout << "method2: = " << dist << ", point_index = " << point_index/D << std::endl;
        dist *= mul_factor;

        if (dist < res.front().first) {
          std::pop_heap (res.begin(), res.end());
          res.pop_back();
          res.push_back({dist, leaves[leaf_index][j]});
          std::push_heap (res.begin(), res.end());
        }
      }
    } else {
      // node
      FT val_q = q[nodes[index].first];
      FT val_n = nodes[index].second;

      // Which child branch should we take first?
      FT diff = val_q - val_n;
      size_t best_child_i, other_child_i;
      if(diff < 0) {
        best_child_i = 2 * index + 1;
        other_child_i = 2 * index + 2;
      }
      else {
        best_child_i = 2 * index + 2;
        other_child_i = 2 * index + 1;
      }

      /*
       *  Partial Euclidean distance, using just one dimension. This is used by the
       *  RKD-tree when computing partial distances while traversing the tree.
       *
       *  Square root is omitted for efficiency.
       */
      float new_dist = mindist + ((val_q - val_n) * (val_q - val_n));

      // collect branch not taken for future use
      if(new_dist * mul_factor < res.front().first) {
        branch.push(std::make_tuple(new_dist, other_child_i, tree_i));
      }

      nn(best_child_i, q, mul_factor, res, c,
         max_leaf_check, vis, branch, mindist, tree_i,
         squared_coords, q_squared_coords, q_index);
    }
  }

  /**
   * \brief Search for nearest neighbors.
   *
   * @param q                - original query
   * @param mul_factor       - approximation factor
   * @param max_leaf_check   - max number of
   * leaves to check
   * @param res              - vector of results (NNs)
   * @param branch           - minimum heap of
   * unvisited branches
   * @param mindist          - current minimum distance
   * @param node_i           - node index
   * @param tree_i           - tree index
   * @param visited          - points visited
   * @param squared_coords   - squared coordinates of
   * points
   * @param q_squared_coords - squared coordinates of
   * queries
   * @param q_index          - index of query
   */
  void search_nn(const std::vector<FT>& q, const float& mul_factor,
                 const int max_leaf_check,
                 std::vector <std::pair<float, int> >& res,
                 Min_heap& branch, const float& mindist,
                 const size_t node_i, const size_t& tree_i,
                 size_t* visited,
                 const std::vector<FT>& squared_coords,
                 const std::vector<FT>& q_squared_coords, int q_index) {
    // edw benw
    int c = 0; // number of leaves visited
    if (rot.size()) {
      std::vector<FT> q_rot;
      householder_transform(1, D, q_rot, q, rot);
      nn(node_i, q_rot, q, mul_factor, res, c, max_leaf_check, visited,
         branch, mindist, tree_i, squared_coords, q_squared_coords, q_index);
    } else {
      nn(node_i, q, mul_factor, res, c, max_leaf_check, visited,
         branch, mindist, tree_i, squared_coords, q_squared_coords, q_index);
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  /////////////////////////// Radius Search /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  /**
   * \brief Recursively perform radius
   * search for nearest neighbors.
   *
   * @param index           - index of node/leaf we are accessing
   * @param q_rot           - rotated query
   * @param q               - original query
   * @param r_epsilon       - approximate radius
   * @param res             - vector of results (NNs)
   * @param c               - number of leaves visited
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param r               - radius
   * @param vis             - points already visited
   */
  void nn(size_t index, const std::vector<FT>& q_rot, const std::vector<FT>& q,
          const float& r_epsilon, std::vector <std::pair<int, float> >& res,
          int& c, const int max_leaf_check, const float& r,
          std::unordered_set<size_t>& vis) {
    if(c ==  max_leaf_check)
      return;
    if (index >= n) {
      ++c;
      // leaf
      size_t point_index, leaf_index = index - n;
      float dist;
      for (size_t j = 0; j < nmax; ++j) {
        // leaves[i][j] gives you the index of the j-th point,
        // located in the i-th leaf.
        // (*points)[leaves[leaf_index][j]] is where the first
        // coordinate of the point described above lies.
        point_index = leaves[leaf_index][j] * D;
        if(vis.find(point_index) != vis.end()) {
          continue;
        }
        vis.insert(point_index);
        dist = Eucl_distance(*points, point_index, point_index + D, q);

        if (dist < r) {
          res.push_back({leaves[leaf_index][j], dist});
        }
      }
    } else {
      // node
      if (q_rot[nodes[index].first] < nodes[index].second) {
        // go to left child
        nn(2 * index + 1, q_rot, q, r_epsilon, res, c, max_leaf_check, r, vis);
        if (q_rot[nodes[index].first] + r_epsilon >= nodes[index].second) {
          nn(2 * index + 2, q_rot, q, r_epsilon, res, c, max_leaf_check, r, vis);
        }
      } else {
        // go to right child
        nn(2 * index + 2, q_rot, q, r_epsilon, res, c, max_leaf_check, r, vis);
        if (q_rot[nodes[index].first] - r_epsilon < nodes[index].second) {
          nn(2 * index + 1, q_rot, q, r_epsilon, res, c, max_leaf_check, r, vis);
        }
      }
    }
  }

  /**
   * \brief Recursively perform radius
   * search for nearest neighbors.
   *
   * @param index           - index of node/leaf we are accessing
   * @param q               - original query
   * @param r_epsilon       - approximate radius
   * @param res             - vector of results (NNs)
   * @param c               - number of leaves visited
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param r               - radius
   * @param vis             - points already visited
   */
  void nn(size_t index, const std::vector<FT>& q, const float& r_epsilon,
          std::vector <std::pair<int, float> >& res, int& c,
          const int max_leaf_check, const float& r,
          std::unordered_set<size_t>& vis) {
    if(c == max_leaf_check)
      return;
    if (index >= n) {
      ++c;
      // leaf
      size_t point_index, leaf_index = index - n;
      float dist;
      for (size_t j = 0; j < nmax; ++j) {
        // leaves[i][j] gives you the index of the j-th point,
        // located in the i-th leaf.
        // (*points)[leaves[leaf_index][j]] is where the first
        // coordinate of the point described above lies.
        point_index = leaves[leaf_index][j] * D;
        if(vis.find(point_index) != vis.end()) {
          continue;
        }
        vis.insert(point_index);
        dist = Eucl_distance(*points, point_index, point_index + D, q);

        if (dist < r) {
          res.push_back({leaves[leaf_index][j], dist});
        }
      }
    } else {
      // node
      if (q[nodes[index].first] < nodes[index].second) {
        // go to left child
        nn(2 * index + 1, q, r_epsilon, res, c, max_leaf_check, r, vis);
        if (q[nodes[index].first] + r_epsilon >= nodes[index].second) {
          nn(2 * index + 2, q, r_epsilon, res, c, max_leaf_check, r, vis);
        }
      } else {
        // go to right child
        nn(2 * index + 2, q, r_epsilon, res, c, max_leaf_check, r, vis);
        if (q[nodes[index].first] - r_epsilon < nodes[index].second) {
          nn(2 * index + 1, q, r_epsilon, res, c, max_leaf_check, r, vis);
        }
      }
    }
  }

  /**
   * \brief Radius search for nearest neighbors.
   *
   * @param q               - original query
   * @param mul_factor      - approximation factor
   * @param max_leaf_check  - max number of
   * leaves to check
   * @param r               - radius
   * @param res             - vector of results (NNs)
   * @param visited         - visited points
   */
  void search_nn_radius(const std::vector<FT>& q, float& mul_factor, int max_leaf_check,
                        float r, std::vector <std::pair<int, float> >& res,
                        std::unordered_set<size_t>& visited) {
    int c = 0; // number of leaves visited
    float r_epsilon = r * mul_factor;
    if (rot.size()) {
      std::vector<FT> q_rot;
      householder_transform(1, D, q_rot, q, rot);
      nn(0, q_rot, q, r_epsilon, res, c, max_leaf_check, r, visited);
    } else {
      nn(0, q, r_epsilon, res, c, max_leaf_check, r, visited);
    }
  }

  /**
   * \brief Get the points of the tree.
   *
   * @return - the points
   */
  const std::vector<FT>* get_points() {
    return points;
  }

  /**
   * \brief Get the dimension of
   * the points of the tree.
   *
   * @return - the dimension
   */
  const size_t dim() {
    return D;
  }

 private:
  /**
   * Coordinates of points
   */
  const std::vector<FT>* points;
  /**
   * Dimension of points
   */
  size_t D;
  /**
   * Leaves of tree
   */
  std::vector<std::vector<size_t> > leaves;  // [l][nmax]
  /**
   * Number of leaves
   */
  size_t l;
  /**
   * Max number of points per leaf
   */
  size_t nmax;
  /**
   * Nodes of tree
   */
  std::vector<std::pair<size_t, FT> > nodes;  // [n]
  /**
   * Number of nodes
   */
  size_t n;
  /**
   * Rotation vector
   */
  std::vector<float> rot;     // rotation vector
};

#endif
