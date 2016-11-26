/** \example main_no_par.cpp
 * This is an example of how to use serial building of the forest
 * and then search it efficiently.
 */
 
/**
@mainpage

This project provides approximate and exact nearest neighbor search for high dimensions.

@author Georgios Samaras
@date 22/04/2015
@version 1.0
*/

#include <string>
#include "../source/Auto_random_kd_forest.h"

int main(int argc, char *argv[]) {
  size_t N = 11, D = 10, Q = 1;
  int k = 2;
  double epsilon = 0.0;
  std::string datafile = "test_files/data.txt", queryfile = "test_files/query.txt";
  std::vector< std::vector<std::pair<float, int> > > res;

  Params mypars;
  mypars.points_per_leaf = 1;
  mypars.trees_no = 1;
  mypars.t = 1;
  mypars.max_leaf_check = 3;
  mypars.rotate_option = No;
  mypars.shuffle_enable = false;
  mypars.sample_size = N;

  Auto_random_kd_forest<float> RKDf(N, D, datafile, Q, queryfile, k, epsilon,
                                    res, &mypars);

  std::cout << "\nRESULTS\n";
  for (std::vector<std::pair<float, int> >::const_iterator it = res[0]
      .begin(); it != res[0].end(); ++it)
    std::cout << it->first << ", index = " << it->second << std::endl;

  std::cout << "main_no_par successfully terminated" << std::endl;
  return 0;
}
