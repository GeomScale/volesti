/** \example main_par_auto.cpp
 * This is an example of how to use parallel building of the forest
 * and then search it efficiently, by using auto configured parameters.
 */

#include <string>
#include "../source/Auto_random_kd_forest.h"

int main(int argc, char *argv[]) {
  size_t N = 11, D = 10, Q = 1;
  int k = 2;
  double epsilon = 0.0;
  std::string datafile = "test_files/data.txt", queryfile = "test_files/query.txt";
  std::vector< std::vector<std::pair<float, int> > > res;

  Auto_random_kd_forest<int> RKDf(N, D, datafile, Q, queryfile, k, epsilon, res);
	
  std::cout << "\nRESULTS\n";
  for (std::vector<std::pair<float, int> >::const_iterator it = res[0]
      .begin(); it != res[0].end(); ++it)
    std::cout << it->first << ", index = " << it->second << std::endl;

		
  std::cout << "main_par_auto successfully terminated" << std::endl;
  return 0;
}
